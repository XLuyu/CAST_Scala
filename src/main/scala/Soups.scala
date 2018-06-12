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

  def support(a:Array[Array[Double]],b:Array[Array[Double]]) = {
    (for ( i <- a.indices ) yield{
      val q = a(i).indices.filter(_!=i).map(j=>a(i)(j))
      val p = b(i).indices.filter(_!=i).map(j=>b(i)(j))
      Math.abs(Util.PearsonCorrelationSimilarity(q,p))
    }).sum/a.length
  }
  def getOneChromGVFromScanner:ArrayBuffer[(String,(Matrix,Matrix))] ={
    // scan this contig to get all snp sites
    val contig = bamScanner.getChrom
    val contiglen = bamScanner.header.getSequence(contig).getSequenceLength
    var sites = new ArrayBuffer[(Int,GenotypeVector)](contiglen/1000)
    do {
      val gv = new GenotypeVector(bamScanner.get())
      if ( gv.isReliable && gv.isHeterogeneous) {
        sites.append((bamScanner.pos, gv))
      }
    } while (bamScanner.nextPosition())
    val coverage = bamScanner.getCoverageLimit
    coverage.foreach(x=>print(f"$x%.0f, "))
    sites = sites.filter { case (_,gv) => gv.consistentWithDepth(coverage) }
    if (sites.length<2) {
      println(f"\n[$contig] SNP: ${sites.length} sites (No enough SNP)")
      return new ArrayBuffer[(String,(Matrix,Matrix))]()
    } else {
      println(f"\n[$contig] SNP: ${sites.length} sites")
    }
    // train matrix
    var headMatrix,tailMatrix = new Matrix(size=bamFilenames.length)
    val snapshot = new ArrayBuffer[Matrix](sites.length)
    val segment = new ArrayBuffer[(String,(Matrix,Matrix))]()
    var lastpos = 0
    for ( (pos,gv) <- sites){
      val matrix = new Matrix(gv.toDistMatrix)
      val decayN = Math.pow(decay,pos-lastpos)
      tailMatrix *= decayN += (matrix *= ((1-decayN)/(1-decay)))
      snapshot.append(tailMatrix.copy)
      lastpos = pos
    }
    lastpos = contiglen + 1
    var lastSegEnd = sites.length-1
    for ( i <- sites.indices.reverse){
      val (pos,gv) = sites(i)
      val matrix = new Matrix(gv.toDistMatrix)
      val decayN = Math.pow(decay,lastpos-pos)
      headMatrix *= decayN += (matrix *= ((1-decayN)/(1-decay)))
      if (0<i && support(snapshot(i-1).normalized.data,headMatrix.copy.normalized.data)<0.3){
        val seg = (f"$contig[${sites(i-1)._1}~${sites(i)._1}:"+
          f"${if (lastSegEnd==sites.length-1) contiglen else sites(lastSegEnd)._1}~"+
          f"${if (lastSegEnd==sites.length-1) contiglen else sites(lastSegEnd+1)._1}]",
          (headMatrix.copy.normalized,snapshot(lastSegEnd).copy.normalized))
        segment.append(seg)
        lastSegEnd = i-1
      }
    }
    val seg = (f"$contig[1~1:"+
      f"${if (lastSegEnd==sites.length-1) contiglen else sites(lastSegEnd)._1}~"+
      f"${if (lastSegEnd==sites.length-1) contiglen else sites(lastSegEnd+1)._1}]",
      (headMatrix.copy.normalized,snapshot(lastSegEnd).copy.normalized))
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
