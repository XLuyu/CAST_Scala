import java.io.File
import java.util
import org.tc33.jheatchart.HeatChart
import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer
import htsjdk._
import htsjdk.samtools.SAMRecord
import me.tongfei.progressbar._
import org.jgrapht.alg.matching.MaximumWeightBipartiteMatching
import org.jgrapht.graph._

class Soups(val bamFilenames:Array[String]){
  val bamScanner = new BamFilesParallelScanner(bamFilenames)
  val stringer = new Stringer()
  def checkAndNormalize(v:Array[Int]): Array[Double] ={
    val sum = v.sum.toDouble
//    val v2 = v.map(x=>if (x>=sum*0.9) 1 else 0 )
//    if (v2.sum!=1 || sum<5) null else v2.map(_.toDouble)
    if (sum<5) null else v.map(_/sum)
  }
  def Normalize(v:Array[Double]): Array[Double] ={
    val sum = v.max
    if (sum==0) v else v.map(_/sum)
  }
  def PearsonCorrelationSimilarity(vec1 : Vector[Double], vec2 : Vector[Double]): Double = {
    val sum_vec1 = vec1.sum
    val sum_vec2 = vec2.sum

    val square_sum_vec1 = vec1.map(x => x * x).sum
    val square_sum_vec2 = vec2.map(x => x * x).sum

    val zipVec = vec1.zip(vec2)

    val product = zipVec.map(x => x._1 * x._2).sum
    val numerator = product - (sum_vec1 * sum_vec2 / vec1.length)

    val dominator = math.sqrt((square_sum_vec1 - sum_vec1*sum_vec1 / vec1.length) * (square_sum_vec2 - sum_vec2*sum_vec2 / vec2.length))

    if(dominator == 0)
      0
    else
      numerator / (dominator * 1.0)
  }

  def genotypeVectorToDistanceMatrix(v:Array[Genotype]) = {
    for ( gi <- v) yield
      for ( gj <- v) yield {
        val d = gi.distance(gj)
        if (d < 0.05) 0 else d
      }
  }
  def getOneChromGVFromScanner:(String,(Matrix,Matrix)) ={
    // scan this contig to get all snp sites
    val contig = bamScanner.getChrom
    val contiglen = bamScanner.header.getSequence(contig).getSequenceLength
    val sites = ArrayBuffer[(Int,Array[Genotype])]()
    do {
      val genotypeVector = bamScanner.get()
      val distMatrix = genotypeVectorToDistanceMatrix(genotypeVector)
      if (!genotypeVector.exists(_.isUnreliable) && distMatrix.flatten.exists(_>=0.05)) {
        sites.append((bamScanner.pos, genotypeVector))
//        println(bamScanner.pos)
      }
    } while (bamScanner.nextPosition())
    if (sites.length<2) {
      println(f"No enough SNP, only ${sites.length} site")
      return (contig,null)
    }
    val intervals = new ArrayBuffer[(Int,Int,Array[Genotype])](contiglen/1000)
    intervals.append((sites(0)._1,(sites(0)._1+sites(1)._1)/2,sites(0)._2))
    for ( i <- 1 until sites.length-1) intervals.append((
      (sites(i-1)._1+sites(i)._1)/2+1,
      (sites(i+1)._1+sites(i)._1)/2,
      sites(i)._2))
    intervals.append(((sites(sites.length-2)._1+sites.last._1)/2+1,sites.last._1,sites.last._2))
    // train matrix
    var headMatrix,tailMatrix = new Matrix(bamFilenames.length)
    for ( (start,end,vector) <- intervals){
      val matrix = genotypeVectorToDistanceMatrix(vector)
      for ( i <- start to end) tailMatrix *= 0.99995 += matrix
    }
    for ( (start,end,vector) <- intervals.reverseIterator){
      val matrix = genotypeVectorToDistanceMatrix(vector)
      for ( i <- end to start by -1) headMatrix *= 0.99995 += matrix
    }
    headMatrix.normalized
    tailMatrix.normalized
    new HeatChart(headMatrix.data).saveToFile(new File(s"png/$contig L.png"))
    new HeatChart(tailMatrix.data).saveToFile(new File(s"png/$contig R.png"))
    (contig,(headMatrix,tailMatrix))
  }
  def getAllGVFromScanner ={
    val GVmap = mutable.Map[String,(Matrix,Matrix)]()
    while (bamScanner.hasNext) {
      val CGV = getOneChromGVFromScanner
      if (CGV._2!=null) GVmap(CGV._1) = CGV._2
    }
    GVmap
  }
  def pairwiseMutualBest(contigsGV: mutable.Map[String, Matrix]) = {
    def support(a:Array[Array[Double]],b:Array[Array[Double]]) = {
//      var result = 0.0
//      for ( i <- a.indices ; j <- a(i).indices)
//        result += (a(i)(j) - b(i)(j))*(a(i)(j) - b(i)(j))
//      result
      (for ( i <- a.indices ) yield{
        val q = a(i).indices.filter(_!=i).map(j=>a(i)(j))
        val p = b(i).indices.filter(_!=i).map(j=>b(i)(j))
        PearsonCorrelationSimilarity(q.toVector,p.toVector)
      }).sum/a.length
    }

    val bestMatch = for ( (k,v) <- contigsGV ) yield {
      val contig = k.substring(1)
//      val similarities = contigsGV.filter(_._1.substring(1)!=contig).mapValues(v2=>support(v.data,v2.data))
      val similarities = contigsGV.filter(_._1!=k).mapValues(v2=>support(v.data,v2.data))
      val best = similarities.maxBy(_._2)
      if (similarities.count(_._2==best._2)>1) (k,null) else (k,best)
    }
//    bestMatch.toArray.sortBy(_._1).foreach(println)
    bestMatch.filter { case (k,v) => v!=null && bestMatch(v._1)._1==k }
  }
  def run(): Unit ={
//    val coverages = bamScanner.getCoverage
//    println(coverages.mkString(" "))
    val contigsTwoEndGV = getAllGVFromScanner
    val contigsOneEndGV = contigsTwoEndGV.flatMap {case (contig,(v1,v2)) => List(("+"+contig,v1),("-"+contig,v2))}
    val finalLink = pairwiseMutualBest(contigsOneEndGV)
    for ( (k,v) <- finalLink if k<v._1) {
      println(k,v)
      stringer.put(k,v._1,v._2)
    }
    stringer.getPath.foreach(println)
  }
}

object Soups {
  def main(args:Array[String]): Unit = new Soups(args).run()
}
