import htsjdk.samtools
import htsjdk.samtools.SAMRecord
import me.tongfei.progressbar.ProgressBar
import scala.collection.JavaConverters._
import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable

class BamFileScanner(filename:String){
  val bamFile = samtools.SamReaderFactory.makeDefault().open(new java.io.File(filename))
  val header = bamFile.getFileHeader.getSequenceDictionary
  val fileIterator = bamFile.iterator.asScala.buffered
  var spanReads = ArrayBuffer[SAMRecord]()
  var sw = new SlidingWindowArray(100000)
  var badRegion = new SlidingWindowInt(100000)

  val code = Map[Char,Int]('A'->0,'C'->1,'G'->2,'T'->3)
  var maxReadSpan, cached = 0
  var cigarRegex = raw"(\d+)[SH]".r
  var coverageStat = mutable.Map[Int,Long]().withDefaultValue(0)

  def isBadRead(read:SAMRecord) = read.getMappingQuality<30 || read.getMateReferenceName!=read.getReferenceName ||
                                  cigarRegex.findAllMatchIn(read.getCigarString).map(_.group(1).toInt).sum+read.getIntegerAttribute("NM")>0.1*read.getReadLength ||
                                  read.hasAttribute("XA")

  private def updateSpanReadsByPosition(cid:Int,pos:Int): Unit ={
    spanReads = spanReads.filter(record=>record.getReferenceIndex==cid && pos<=record.getEnd)
    cached = pos + sw.size/2
    var continue = true
    while (fileIterator.hasNext && continue )
      if (fileIterator.head.getReadUnmappedFlag)
        fileIterator.next()
      else {
        val icid = fileIterator.head.getReferenceIndex
        val ipos = fileIterator.head.getAlignmentStart
        if (icid<cid) fileIterator.next() else
        if (icid==cid && ipos<=cached){
          val record = fileIterator.next()
          val bad = isBadRead(record)
          maxReadSpan = Math.max(maxReadSpan,record.getLengthOnReference)
          if (record.getLengthOnReference>sw.size/3)
            spanReads.append(record)
          else {
            val read = record.getReadString
            for ( i <- 1 to record.getReadLength if code.contains(read(i-1)) && record.getReferencePositionAtReadPosition(i)>0){
              sw.inc(record.getReferencePositionAtReadPosition(i),code(read(i-1)))
              if (bad) badRegion.inc(record.getReferencePositionAtReadPosition(i))
            }
          }
        }
        else continue = false
      }
  }
  def get(cid:Int,pos:Int) = {
    if (pos==1) {
      sw = new SlidingWindowArray(100000)
      badRegion = new SlidingWindowInt(100000)
      cached = 0
    }
    if (pos>cached) updateSpanReadsByPosition(cid,pos)
    val vector = Array(0,0,0,0)
    var badCount = 0
    for ( record <- spanReads) {
      val rpos = record.getReadPositionAtReferencePosition(pos)
      if (rpos>0) {
        val c = record.getReadString()(rpos-1)
        if (code.contains(c)) vector(code(c)) += 1
        if (isBadRead(record)) badCount += 1
      }
    }
    badCount += badRegion.getAndClean(pos)
    val count = (vector zip sw.getAndClean(pos)).map(x=>x._1+x._2)
    val genotype = new Genotype(count.map(_.toDouble),badCount)
    if (genotype.isNonRepeat) coverageStat(count.sum) += 1
    genotype
  }
  def getCoverageLimit = {
//    println(f"${coverageStat.foldLeft(0L)((a,b)=>a+b._1*b._2)}/${coverageStat.foldLeft(0)(_+_._2)}")
    val lowerbound = coverageStat.foldLeft(0L)((a,b)=>a+b._1*b._2)/coverageStat.foldLeft(0L)(_+_._2)
    if (coverageStat.isEmpty) 0.0 else coverageStat.filter(_._1>=lowerbound).maxBy(_._2)._1.toDouble
  }
  def getTotalMappedBase ={
    val iter = samtools.SamReaderFactory.makeDefault().open(new java.io.File(filename)).iterator().asScala.buffered
    (for ( i <- iter if !i.getReadUnmappedFlag) yield i.getReadLength.toLong).sum
  }
}
class BamFilesParallelScanner(filenames:Array[String]){
  val bamScanners = filenames.map{filename => new BamFileScanner(filename)}
  val bamScannersPar = bamScanners.par
  val pool = new collection.parallel.ForkJoinTaskSupport(new java.util.concurrent.ForkJoinPool())
  bamScannersPar.tasksupport = pool
  val headers = bamScanners.map(_.header)
  val header = headers(0)
  if (headers.exists(_!=header)) throw new Exception("[Error] Headers in input files are inconsistent")
  def maxReadSpan = bamScanners.map(_.maxReadSpan).max
  var cid = 0
  private var maxpos = header.getSequence(cid).getSequenceLength
  var pos = 0
  val pb = new ProgressBar("Scan Bam files", header.getReferenceLength)
  pb.start()
  nextPosition()

  def nextPosition() = { // return true if on same contig
    pos += 1
    if (pos>maxpos){
      pb.stepBy(pos%1000).setExtraMessage(s"finished $getChrom")
      cid += 1
      pos = 1
      if (cid<header.size) maxpos = header.getSequence(cid).getSequenceLength
      else pb.stop()
    }
    if (pos%1000==0) pb.stepBy(1000)
    pos!=1
  }
  def getCoverage = {
    bamScanners.map{_.getTotalMappedBase/header.getReferenceLength.toDouble}
  }
  def get() = {
    if (pos==1 || pos>bamScanners.head.cached)
      bamScannersPar.map { _.get(cid, pos)}.toArray
    else
      bamScanners.map { _.get(cid, pos)}
  }
  def getCoverageLimit = bamScanners.map(_.getCoverageLimit)
  def getChrom = header.getSequence(cid).getSequenceName
  def hasNext = cid<header.size
}
