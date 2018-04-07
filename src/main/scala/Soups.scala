import java.io.File
import java.util

import scala.collection.JavaConverters._
import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer
import htsjdk._
import htsjdk.samtools.SAMRecord
import me.tongfei.progressbar._
import org.jgrapht._
import org.jgrapht.alg.matching.MaximumWeightBipartiteMatching
import org.jgrapht.graph._

import scala.reflect.ClassTag

class BufferedIterator[A: ClassTag](iter:Iterator[A], size:Int=1000) extends Iterator[A]{
  val buffer = new Array[A](size)
  var tail = 0
  var i = 0

  def head = {
    if (i==tail) reload()
    if (i==tail) throw new Exception("Empty Iterator!")
    buffer(i)
  }

  private def reload(): Unit ={
    i = 0
    tail = 0
    while (iter.hasNext && tail!=buffer.size){
      buffer(tail) = iter.next()
      tail += 1
    }
  }
  override def hasNext = iter.hasNext
  override def next() = {
    if (i==tail) reload()
    if (i==tail) throw new Exception("Empty Iterator!")
    i += 1
    buffer(i-1)
  }
}
class SlidingWindow(val size:Int){
  val window = Array.fill(size)(Array(0,0,0,0))
  def inc(i:Int,j:Int) = window(i%size)(j) += 1
  def getAndClean(i:Int) ={
    val s = window(i%size)
    window(i%size) = Array(0,0,0,0)
    s
  }
}
class BamFileScanner(filename:String){
  val bamFile = samtools.SamReaderFactory.makeDefault().open(new java.io.File(filename))
  val header = bamFile.getFileHeader.getSequenceDictionary
  val fileIterator = bamFile.iterator.asScala.buffered
  var spanReads = ArrayBuffer[SAMRecord]()
  var sw = new SlidingWindow(10000)
  var maxReadSpan = 0

  val code = Map[Char,Int]('A'->0,'C'->1,'G'->2,'T'->3)
  var cid = 0
  private var maxpos = header.getSequence(cid).getSequenceLength
  var pos = 0

  private def updateSpanReadsByPosition(): Unit ={
    spanReads = spanReads.filter(record=>record.getReferenceIndex==cid && pos<=record.getEnd)
    var continue = true
    while (fileIterator.hasNext && continue ){
      if (fileIterator.head.getReadUnmappedFlag || fileIterator.head.getMappingQuality==0) fileIterator.next()
      else {
        val icid = fileIterator.head.getReferenceIndex
        val ipos = fileIterator.head.getAlignmentStart
        if (icid<cid || icid==cid && ipos<=pos){
          val record = fileIterator.next()
          maxReadSpan = Math.max(maxReadSpan,record.getLengthOnReference)
          if (record.getLengthOnReference>sw.size/2)
            spanReads.append(record)
          else {
            val read = record.getReadString
            for ( i <- 1 to record.getReadLength if code.contains(read(i-1)) && record.getReferencePositionAtReadPosition(i)>0){
              sw.inc(record.getReferencePositionAtReadPosition(i),code(read(i-1)))
            }
          }
        }
        else continue = false
      }
    }
  }
  def nextPosition() = { // return true if on same contig
    pos += 1
    if (pos>maxpos){
      cid += 1
      sw = new SlidingWindow(10000)
      if (cid<header.size) maxpos = header.getSequence(cid).getSequenceLength
      pos = 1
    }
    pos!=1
  }
  def get() = {
    updateSpanReadsByPosition()
    val vector = Array(0,0,0,0)
    for ( record <- spanReads) {
      val rpos = record.getReadPositionAtReferencePosition(pos)
      if (rpos>0) {
        val c = record.getReadString()(rpos-1)
        if (code.contains(c)) vector(code(c)) += 1
      }
    }
    (vector zip sw.getAndClean(pos)).map(x=>x._1+x._2)
  }
}
class BamFilesParallelScanner(filenames:Array[String]){
  val bamScanners = filenames.map{filename => new BamFileScanner(filename)}
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
      pb.stepBy(pos%1000)
      pb.setExtraMessage(s"finished $getChrom")
      cid += 1
      if (cid<header.size) maxpos = header.getSequence(cid).getSequenceLength
      else pb.stop()
      pos = 1
    }
    if (pos%1000==0) pb.stepBy(1000)
    bamScanners.map(_.nextPosition())
    pos!=1
  }
  def get() = bamScanners.map { _.get()}
  def getChrom = header.getSequence(cid).getSequenceName
  def hasNext = cid<header.size
}
class BamLoader(){
  val code = Map[Char,Int]('A'->0,'C'->1,'G'->2,'T'->3)
  def contigs = mutable.Map[String,Int]() // get contig length by name
  def loadBaseCountOnGenomeFromBam(bamFilename: String) = {
    val startTimeStamp = System.currentTimeMillis()/1000.0
    val bamFile = htsjdk.samtools.SamReaderFactory.makeDefault().open(new java.io.File(bamFilename))
    val header = bamFile.getFileHeader
    val contigs = header.getSequenceDictionary
    val counter = mutable.Map[String,Array[Array[Int]]]() // counter(contig)(pos)(base) = count
    for ( i <- 0 until contigs.size){
      val contig = contigs.getSequence(i)
      this.contigs(contig.getSequenceName) = contig.getSequenceLength
      counter(contig.getSequenceName) = Array.fill(contig.getSequenceLength+1)(Array(0,0,0,0))
    }
    var recordCount = 0
    val pb = new ProgressBar("Tester", 10000)
    for ( record <- ProgressBar.wrap(bamFile.iterator,s"Loading $bamFilename").asScala) {
      recordCount += 1
      if (recordCount%10000==0) pb.step()
      val read = record.getReadString
      if (!record.getReadUnmappedFlag)
        for ( i <- 1 to record.getReadLength if read(i-1)!='N'){
          counter(record.getContig)(record.getReferencePositionAtReadPosition(i))(code(read(i-1))) += 1
        }
    }
    println(f"[Success] load $recordCount records from $bamFilename in ${System.currentTimeMillis()/1000.0-startTimeStamp}%.1f seconds")
    counter
  }
}

class Soups(val bamFilenames:Array[String]){
  val bamLoader = new BamLoader()
  val bamScanner = new BamFilesParallelScanner(bamFilenames)
  def checkAndNormalize(v:Array[Int]): Array[Double] ={
    val sum = v.sum.toDouble
    val v2 = v.map(x=>if (x>=sum*0.9) 1 else 0 )
    if (sum<5) return null
    if (v2.sum==1) v2.map(_.toDouble) else null
  }
  def distance(a:Array[Double],b:Array[Double]) = if (a==null || b==null) 10.0 else (a zip b map {case (x,y)=>(x-y)*(x-y)}).sum / 2
  def getTwoEndGV(contig: String, length:Int, bamPileUps: Array[mutable.Map[String, Array[Array[Int]]]]):(Array[Int],Array[Int]) = {

    val GenotypesVectorTable:ArrayBuffer[Array[Int]] = ArrayBuffer[Array[Int]]()
    val GenotypesVectorCount:ArrayBuffer[Int] = ArrayBuffer[Int]()
    val GenotypesVectorPosSum:ArrayBuffer[Double] = ArrayBuffer[Double]()
    val contigPileUps = bamPileUps.map(_(contig))
    for ( pos <- 1 to length) {
      val ufs = new UnionFindSet()
      val ratios = contigPileUps.map(_(pos)).map(checkAndNormalize)
      for ( i <- ratios.indices )
        for ( j <- 0 until i if distance(ratios(i),ratios(j))<0.05)
          ufs.union(i,j)
      val GenotypesVector = ratios.indices.map(x=>if (ratios(x)==null) -1 else ufs.find(x)).toArray
      if (GenotypesVector.exists(_!=0)){
        var GVindex = GenotypesVectorTable.indexWhere(_.sameElements(GenotypesVector))
        if (GVindex == -1) {
          GVindex = GenotypesVectorTable.length
          GenotypesVectorTable.append(GenotypesVector)
          GenotypesVectorCount.append(0)
          GenotypesVectorPosSum.append(0.0)
        }
        GenotypesVectorCount(GVindex) += 1
        GenotypesVectorPosSum(GVindex) += pos
      }
    }
    if (GenotypesVectorTable.isEmpty) return null
    if (GenotypesVectorTable.length==1) return (GenotypesVectorTable(0),GenotypesVectorTable(0))
    val mostCommon = GenotypesVectorTable.indices.sortBy(GenotypesVectorCount)
    val top = mostCommon(0)
    val second = mostCommon(1)
    if (GenotypesVectorPosSum(top)/GenotypesVectorCount(top)<GenotypesVectorPosSum(second)/GenotypesVectorCount(second))
      (GenotypesVectorTable(top),GenotypesVectorTable(second))
    else
      (GenotypesVectorTable(second),GenotypesVectorTable(top))
  }

  def pairwiseMutualBest(contigsGV: mutable.Map[String, Array[Int]]) = {
    def support(x:Array[Int],y:Array[Int]) = {
      val graph = new SimpleWeightedGraph[Int,DefaultWeightedEdge](classOf[DefaultWeightedEdge])
      val X = new java.util.HashSet[Int]()
      val Y = new java.util.HashSet[Int]()
      val edges = (x zip y) filter { case (a,b) => a>=0 && b>=0 } groupBy identity mapValues(_.length)
      for ( ((s,t),w) <- edges ) {
        graph.addVertex(s)
        X.add(s)
        graph.addVertex(-t-1)
        Y.add(-t-1)
        graph.setEdgeWeight(graph.addEdge(s,-t-1),w)
      }
      val mwbm = new MaximumWeightBipartiteMatching(graph,X,Y)
      mwbm.getMatching.getWeight.round.toInt
    }

    val bestMatch = for ( (k,v) <- contigsGV ) yield {
      val contig = k.substring(1)
      val similarities = contigsGV.filter(_._1.substring(1)!=contig).mapValues(v2=>support(v,v2))
      val best = similarities.maxBy(_._2)
      if (similarities.count(_._2==best._2)>1) (k,"") else (k,best._1)
    }
    bestMatch.filter { case (k,v) => v!="" && bestMatch(v)==k }
  }
  def getAllGVFromScanner ={
    val GVmap = mutable.Map[String,(Array[Int],Array[Int])]()
    while (bamScanner.hasNext) {
      val CGV = getOneChromGVFromScanner
      if (CGV._2!=null) GVmap(CGV._1) = CGV._2
    }
    GVmap
  }
  def getOneChromGVFromScanner:(String,(Array[Int],Array[Int])) ={
    class GVSpan(val GV:Array[Int], pos:Int, var count:Int=1){
      var left,right = pos
      override def toString = s"${GV.mkString(",")}\t[$left--($count)-->$right]"
    }
    val GVSpanTable:ArrayBuffer[GVSpan] = ArrayBuffer[GVSpan]()
    val contig = bamScanner.getChrom
    do {
      val ufs = new UnionFindSet()
      val ratios = bamScanner.get().map(checkAndNormalize).toArray
      for ( i <- ratios.indices )
        for ( j <- 0 until i if distance(ratios(i),ratios(j))<0.05)
          ufs.union(i,j)
      val GV = ratios.indices.map(x=>if (ratios(x)==null) -1 else ufs.find(x)).toArray
      if (!GV.contains(-1) && GV.filter(_>=0).distinct.length>1){
        if (GVSpanTable.nonEmpty &&
          !GVSpanTable.last.GV.sameElements(GV) &&
          GVSpanTable.last.right-GVSpanTable.last.left<bamScanner.maxReadSpan){
          GVSpanTable.trimEnd(1)
        }
        if (GVSpanTable.nonEmpty && GVSpanTable.last.GV.sameElements(GV)){
          GVSpanTable.last.count += 1
          GVSpanTable.last.right = bamScanner.pos
        } else {
          GVSpanTable.append(new GVSpan(GV,bamScanner.pos))
        }
      }
    } while (bamScanner.nextPosition())
    if (GVSpanTable.nonEmpty &&
      GVSpanTable.last.right-GVSpanTable.last.left<bamScanner.maxReadSpan){
      GVSpanTable.trimEnd(1)
    }
    if (GVSpanTable.isEmpty) return (contig,null)
    if (GVSpanTable.length==1) return (contig,(GVSpanTable(0).GV,GVSpanTable(0).GV))
    println()
    println(contig,"\t",GVSpanTable.head.toString())
    println(contig,"\t",GVSpanTable.last.toString())
    if (bamScanner.cid==1)
      for (i <- GVSpanTable)
        println(i.toString)
    (contig,(GVSpanTable.head.GV,GVSpanTable.last.GV))
  }
  def run(): Unit ={
    val contigsTwoEndGV = getAllGVFromScanner
    val contigsOneEndGV = contigsTwoEndGV.flatMap {case (contig,(v1,v2)) => List(("+"+contig,v1),("-"+contig,v2))}
    val finalLink = pairwiseMutualBest(contigsOneEndGV)
    for ( (k,v) <- finalLink if k<v) println(k,v)
  }
}

object Soups {
  def main(args:Array[String]): Unit = new Soups(args).run()
}
