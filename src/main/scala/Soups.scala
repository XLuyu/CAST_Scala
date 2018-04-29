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
  val decay = 0.99995
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
  def PearsonCorrelationSimilarity(vec1 : IndexedSeq[Double], vec2 : IndexedSeq[Double]): Double = {
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
  def support(a:Array[Array[Double]],b:Array[Array[Double]]) = {
    //      var result = 0.0
    //      for ( i <- a.indices ; j <- a(i).indices)
    //        result += (a(i)(j) - b(i)(j))*(a(i)(j) - b(i)(j))
    //      result
    (for ( i <- a.indices ) yield{
      val q = a(i).indices.filter(_!=i).map(j=>a(i)(j))
      val p = b(i).indices.filter(_!=i).map(j=>b(i)(j))
      PearsonCorrelationSimilarity(q,p)
    }).sum/a.length
  }
  def getOneChromGVFromScanner:ArrayBuffer[(String,(Matrix,Matrix))] ={
    // scan this contig to get all snp sites
    val contig = bamScanner.getChrom
    val contiglen = bamScanner.header.getSequence(contig).getSequenceLength
    val maxMisassembleSpan = bamScanner.header.getReferenceLength/100
    val sites = new ArrayBuffer[(Int,Array[Genotype])](contiglen/1000)
    do {
      val genotypeVector = bamScanner.get()
      if ( !genotypeVector.exists(_.isUnreliable) && //contig=="tig00000046_pilon" &&
        !genotypeVectorToDistanceMatrix(genotypeVector).exists(x=>{ val sx = x.sorted; sx.last-sx(1)<0.2})) {
        sites.append((bamScanner.pos, genotypeVector))
      }
    } while (bamScanner.nextPosition())
    if (sites.length<2) {
      println(f"\n[$contig] SNP: ${sites.length} sites (No enough SNP)")
      return new ArrayBuffer[(String,(Matrix,Matrix))]()
    } else {
      println(f"\n[$contig] SNP: ${sites.length} sites")
    }
    val intervals = new ArrayBuffer[(Int,Int,Array[Genotype])](sites.length)
    intervals.append((sites(0)._1,(sites(0)._1+sites(1)._1)/2,sites(0)._2))
    for ( i <- 1 until sites.length-1) intervals.append((
      (sites(i-1)._1+sites(i)._1)/2+1,
      (sites(i+1)._1+sites(i)._1)/2,
      sites(i)._2))
    intervals.append(((sites(sites.length-2)._1+sites.last._1)/2+1,sites.last._1,sites.last._2))
    // train matrix
    var headMatrix,tailMatrix = new Matrix(bamFilenames.length)
    val snap = new ArrayBuffer[Matrix](intervals.length)
    val segment = new ArrayBuffer[(String,(Matrix,Matrix))]()
    for ( (start,end,vector) <- intervals){
      val matrix = new Matrix(vector.length,genotypeVectorToDistanceMatrix(vector))
      val decayN = Math.pow(decay,end-start+1)
      tailMatrix *= decayN += (matrix *= ((1-decayN)/(1-decay)))
      snap.append(tailMatrix.copy)
//      new HeatChart(tailMatrix.data).saveToFile(new File(f"png/$start%07d~$end%07d R.png"))
    }
    var last,i = intervals.length-1
    for ( (start,end,vector) <- intervals.reverseIterator){
      val matrix = new Matrix(vector.length,genotypeVectorToDistanceMatrix(vector))
      val decayN = Math.pow(decay,end-start+1)
      headMatrix *= decayN += (matrix *= ((1-decayN)/(1-decay)))
      i -= 1
      if (i>0 && //sites(i+1)._1-sites(i)._1<maxMisassembleSpan &&
//        sites(last)._1-sites(i+1)._1>maxMisassembleSpan &&
//        sites(i)._1>maxMisassembleSpan &&
        support(snap(i).normalized.data,headMatrix.copy.normalized.data)<0.3){
        val seg = (f"$contig[${sites(i)._1}~${sites(i+1)._1}:${
          if (last==sites.length-1) contiglen else sites(last)._1}~${
          if (last==sites.length-1) contiglen else sites(last+1)._1}]",
          (headMatrix.copy.normalized,snap(last).copy.normalized))
        segment.append(seg)
        last = i
      }
//      new HeatChart(headMatrix.data).saveToFile(new File(f"png/$start%07d~$end%07d L.png"))
    }
    val seg = (f"$contig[1~1:${
      if (last==sites.length-1) contiglen else sites(last)._1}~${
      if (last==sites.length-1) contiglen else sites(last+1)._1}]",
      (headMatrix.copy.normalized,snap(last).copy.normalized))
    segment.append(seg)
    for ( (contig,(headMatrix,tailMatrix)) <- segment){
      new HeatChart(headMatrix.data).saveToFile(new File(s"png/$contig L.png"))
      new HeatChart(tailMatrix.data).saveToFile(new File(s"png/$contig R.png"))
    }
    segment
  }
  def getAllGVFromScanner ={
    val GVmap = mutable.Map[String,(Matrix,Matrix)]()
    while (bamScanner.hasNext) {
      val CGV = getOneChromGVFromScanner
      for ( (contig,(headMatrix,tailMatrix)) <- CGV)
        GVmap(contig) = (headMatrix,tailMatrix)
    }
    GVmap
  }
  def pairwiseMutualBest(contigsGV: mutable.Map[String, Matrix]) = {
    val candidate = ArrayBuffer[(String,(String,Double))]()
    for ( (k,v) <- contigsGV ) yield {
      val contig = k.substring(1)
      val similarities = contigsGV.filter(_._1.substring(1)!=contig).mapValues(v2=>support(v.data,v2.data))
      similarities.filter(_._2>0.6).foreach(x=>candidate.append((k,x)))
    }
    candidate
  }
  def run(): Unit ={
//    val coverages = bamScanner.getCoverage
//    println(coverages.mkString(" "))
    val contigsTwoEndGV = getAllGVFromScanner
    val contigsOneEndGV = contigsTwoEndGV.flatMap {case (contig,(v1,v2)) => List(("+"+contig,v1),("-"+contig,v2))}
    val finalLink = pairwiseMutualBest(contigsOneEndGV)
    for ( (k,v) <- finalLink if k<v._1) {
      println(k,v)
//      stringer.put(k,v._1,v._2)
    }
//    stringer.getPath.foreach(println)
  }
}

object Soups {
  def main(args:Array[String]): Unit = new Soups(args).run()
}
