import htsjdk.samtools
import htsjdk.samtools.SAMRecord
import me.tongfei.progressbar.ProgressBar

import scala.collection.JavaConverters._
import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer
class BamFileDetector(filename:String) extends BamFileScanner(filename){
  var coverageStat = mutable.Map[Int,Long]().withDefaultValue(0)

  override def get(cid:Int,pos:Int) = {
    if (pos==1) cached = 0
    if (pos>cached) updateSpanReadsByPosition(cid,pos)
    val badCount = badRegion.getAndClean(pos)
    val count = sw.getAndClean(pos)
    val genotype = new Genotype(count.map(_.toDouble),badCount)
    if (genotype.isReliable) coverageStat(count.sum) += 1
    genotype
  }
  def getCoverage = {
    val lowerbound = coverageStat.foldLeft(0.0)((a,b)=>a+b._1*b._2)/coverageStat.foldLeft(0.0)(_+_._2)
    if (coverageStat.isEmpty) 0.0 else coverageStat.filter(_._1>=lowerbound).maxBy(_._2)._1.toDouble
  }
}
class DepthDetector(filenames:Array[String]){
  val bamScanners = filenames.map{filename => new BamFileDetector(filename)}
  val bamScannersPar = bamScanners.par
  val pool = new collection.parallel.ForkJoinTaskSupport(new java.util.concurrent.ForkJoinPool())
  bamScannersPar.tasksupport = pool
  val headers = bamScanners.map(_.header)
  val header = headers(0)
  if (headers.exists(_!=header)) throw new Exception("[Error] Headers in input files are inconsistent")
  val contigNum = header.getSequences.size()
  var cid = (0 until contigNum).maxBy(i=>header.getSequence(i).getSequenceLength)
  private val maxpos = header.getSequence(cid).getSequenceLength
  var pos = 0
  val pb = new ProgressBar("Pre-check", maxpos)

  val getCoverage = {
    pb.start()
    while (pos<maxpos){
      pos += 1
      if (pos==1 || pos>bamScanners.head.cached)
        bamScannersPar.map { _.get(cid, pos)}.toArray
      else
        bamScanners.map { _.get(cid, pos)}
      if (pos%1000==0) pb.stepBy(1000)
    }
    pb.stop()
    bamScanners.map(_.getCoverage)
  }
}
