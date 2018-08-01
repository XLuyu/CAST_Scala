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
  var sw = new SlidingWindowArray(100000)
  var badRegion = new SlidingWindowInt(100000)

  val code = Map[Char,Int]('A'->0,'C'->1,'G'->2,'T'->3)
  var cached = 0
  var cigarRegex = raw"(\d+)[SHDI]".r

  def isBadRead(read:SAMRecord) = read.getMappingQuality<30 || read.getMateUnmappedFlag || read.getMateReferenceName!=read.getReferenceName ||
                                  cigarRegex.findAllMatchIn(read.getCigarString).map(_.group(1).toInt).map(x=>x*x-1).sum+read.getIntegerAttribute("NM")>0.1*read.getReadLength ||
                                  read.hasAttribute("XA")

  def updateSpanReadsByPosition(cid:Int,pos:Int): Unit ={
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
          val read = record.getReadString
          val qual = record.getBaseQualities
          for ( i <- record.getStart to record.getEnd) {
            val pos = record.getReadPositionAtReferencePosition(i)
            if (pos==0) {
              sw.inc(i,4)
            } else {
              val ch = read(pos-1)
              if (code.contains(ch) && qual(pos-1)>=20) sw.inc(i,code(ch))
            }
            if (bad) badRegion.inc(i)
          }
        }
        else continue = false
      }
  }
  def get(cid:Int,pos:Int) = {
    if (pos==1) cached = 0
    if (pos>cached) updateSpanReadsByPosition(cid,pos)
    val badCount = badRegion.getAndClean(pos)
    val count = sw.getAndClean(pos)
    val genotype = new Genotype(count.map(_.toDouble),badCount)
    //if (genotype.isNonRepeat) coverageStat(count.sum) += 1
    genotype
  }
}
class BamFilesParallelScanner(filenames:Array[String]){
  val bamScanners = filenames.map{filename => new BamFileScanner(filename)}
  val bamScannersPar = bamScanners.par
  val pool = new collection.parallel.ForkJoinTaskSupport(new java.util.concurrent.ForkJoinPool())
  bamScannersPar.tasksupport = pool
  val headers = bamScanners.map(_.header)
  val header = headers(0)
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
      pos = 1
      do {
        cid += 1
        if (cid<header.size) {
          maxpos = header.getSequence(cid).getSequenceLength
          if (maxpos<5000) pb.stepBy(maxpos)
        }
        else pb.stop()
      } while (cid<header.size && maxpos<5000)
    }
    if (pos%1000==0) pb.stepBy(1000)
    pos!=1
  }
  def get() = {
    if (pos==1 || pos>bamScanners.head.cached)
      bamScannersPar.map { _.get(cid, pos)}.toArray
    else
      bamScanners.map { _.get(cid, pos)}
  }
  def getChrom = header.getSequence(cid).getSequenceName
  def hasNext = cid<header.size
}
