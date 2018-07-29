import java.io.File

import scala.collection.mutable
import me.tongfei.progressbar.ProgressBar
import htsjdk.samtools.reference._
import org.apache.commons.io.FileUtils

import scala.collection.mutable.ArrayBuffer

class FastaWriter(file:File){
  val writer = new java.io.PrintWriter(file)
  def writeRecord(name:String, seq:String) {
    writer.write(f">$name")
    var cnt = 0 
    for ( c <- seq){
      if (cnt%80==0) writer.write('\n')
      writer.write(c)
      cnt += 1
    }
    writer.write('\n')
  }
}
class Stringer(conf:Conf) {
  final val IDpattern = raw"([^\(]+)\((\d+),(\d+)\)\((\d+),(\d+)\)".r
  val fastaFileName = conf.draftFile()
  val tempRef = new File(conf.blastTmpDir+File.separator+"ref.fasta")
  val tempQry = new File(conf.blastTmpDir+File.separator+"qry.fasta")
  val fasta = new FastaSequenceFile(new File(fastaFileName),false)
  val outFasta = new FastaWriter(conf.outputFasta)
  val contigs = mutable.Map[String,String]()
  class Segment(identifier:String) {
    val groups = IDpattern.findPrefixMatchOf(identifier).get.subgroups.toArray
    var (contig, leftEx, left, right, rightEx) = (groups(0),groups(1).toInt,groups(2).toInt,groups(3).toInt,groups(4).toInt)
    def seq = contigs(contig).substring(leftEx-1,rightEx)
    def len = rightEx-leftEx+1
    var joint = Array[(Link,Int)](null,null)
    def apply(idx:Int) = joint(idx)
    def updateJoint(idx:Int, v:(Link,Int)) {
      joint(idx) = v
      if (v!=null) v._1.joint(v._2) = (this,idx)
    }
  }
  val getSegment = mutable.Map[String,Segment]()
  class OrientedSegment(val segment:Segment, forward:Boolean){
    def this(identifier:String, forward:Boolean=true) = this(getSegment.getOrElseUpdate(identifier,new Segment(identifier)),forward)
    def len = segment.len
    def seq = if (forward) segment.seq else segment.seq.reverseMap(Util.complementary)
    def leftExLen = if (forward) segment.left-segment.leftEx else segment.rightEx-segment.right
    def rightExLen = if (!forward) segment.left-segment.leftEx else segment.rightEx-segment.right
    def pruneL(loss:Int) = if (!forward) segment.rightEx -= loss else segment.leftEx += loss
    def pruneR(loss:Int) = if (forward) segment.rightEx -= loss else segment.leftEx += loss
    def bestL = if (forward) segment(0) else segment(1)
    def bestR = if (forward) segment(1) else segment(0)
    def setBestL(v:(Link,Int)) = if (forward) segment.joint(0)=v else segment.joint(1)=v
    def setBestR(v:(Link,Int)) = if (forward) segment.joint(1)=v else segment.joint(0)=v
    def getName = {
      val l = if (segment.leftEx==1) "S" else segment.leftEx.toString
      val r = if (segment.rightEx==contigs(segment.contig).length) "E" else segment.rightEx.toString
      if (forward) f"${segment.contig}(${l}_$r)" else f"${segment.contig}(${r}_$l)"
    }
  }
  class Contig(){
    var name = new StringBuilder()
    val buffer = new StringBuilder()
    def append(id:String, seq:String){
      name append id
      buffer append seq
    }
  }
  class Link(A:String, B:String, val s:Double){
    val joint = Array((getSegment.getOrElseUpdate(A.substring(1),new Segment(A.substring(1))),if (A(0)=='+') 0 else 1),
                      (getSegment.getOrElseUpdate(B.substring(1),new Segment(B.substring(1))),if (B(0)=='+') 0 else 1))
    def apply(idx:Int, place:Int) = new OrientedSegment(joint(idx)._1,joint(idx)._2!=place)
    var overlap:Array[Int] = null
  }

  val javaRuntime = Runtime.getRuntime
  var records = List[(String,String,Double)]()
  var links:List[Link] = null

  def loadFasta(): Unit ={
    var contig: ReferenceSequence = fasta.nextSequence
    while (contig!=null){
      contigs(contig.getName) = contig.getBaseString.toUpperCase
      contig = fasta.nextSequence
    }
    links = records.map(x=>new Link(x._1,x._2,x._3))
  }
  def BLASTSeqPair(ref:String,qry:String) ={
    FileUtils.write(tempRef,f">reference\n$ref\n")
    FileUtils.write(tempQry,f">query\n$qry\n")
    javaRuntime.exec(f"makeblastdb -in ${tempRef.getAbsolutePath} -dbtype nucl").waitFor()
    val cmd = javaRuntime.exec(Array("blastn", "-db", tempRef.getAbsolutePath, "-query", tempQry.getAbsolutePath, "-outfmt", "6 sstart send qstart qend"))
    scala.io.Source.fromInputStream(cmd.getInputStream).getLines().toArray.map(_.split('\t').map(_.toInt))
    // each line: sstart send qstart qend
  }
  def anchorByBLAST(ref:OrientedSegment, qry:OrientedSegment):Array[Int] = {
    val hsps = BLASTSeqPair(ref.seq,qry.seq)
    val valid = hsps.filter(hsp => hsp(0)<hsp(1) && hsp(2)<hsp(3))
    if (valid.isEmpty) return null
    val best = valid.minBy(x=>Math.max(ref.len-ref.rightExLen-x(1),0)+Math.max(x(2)-1-qry.leftExLen,0))
    val loss = Math.max(ref.len-ref.rightExLen-best(1),0)+Math.max(best(2)-1-qry.leftExLen,0)
    best(0) = ref.len-best(0)+1
    best(1) = ref.len-best(1)+1
    if (loss<=Math.min(ref.len,qry.len)*0.05) best else null
    // best=(dist from ref alignment end to length, dist from ref aligment start to length, qry alignment start, qry alignment end)
  }
  def linkFilter() = {
    val pb = new ProgressBar("BLAST overlap",links.length)
    pb.start()
    for ( link <- links) {
      val ref = link(0,0)
      val qry = link(1,1)
      link.overlap = anchorByBLAST(ref,qry)
      pb.step
      if (link.overlap!=null) {
        if (ref.bestR==null || ref.bestR._1.s<link.s) ref.setBestR(link,0)
        if (qry.bestL==null || qry.bestL._1.s<link.s) qry.setBestL(link,1)
      }
    }
    pb.stop
    println(f"[Info] ${links.count(_.overlap!=null)} candidate links pass loss filter.")
    for ( link <- links if link.overlap!=null) {
      val ref = link(0,0)
      val qry = link(1,1)
      if (ref.bestR==(link,0) && qry.bestL!=(link,1)) ref.setBestR(null)
      if (ref.bestR!=(link,0) && qry.bestL==(link,1)) qry.setBestL(null)
      if (ref.bestR!=(link,0) && qry.bestL!=(link,1)) link.overlap = null
    }
    println(f"[Info] ${links.count(_.overlap!=null)} candidate links pass date filter.")
    val segmentsByContig = getSegment.values.groupBy(_.contig).flatMap { kv=>{
      val segments = kv._2.toArray.sortBy(_.leftEx)
      var last = (1,1)
      var lastJoint:(Link,Int) = null
      var result = List[Segment]()
      for ( segment <- segments) {
        if (segment(0)!=null) {
          if (last != (segment.leftEx, segment.left)) {
            result ::= new Segment(f"${kv._1}$last${(segment.leftEx, segment.left)}")
            result.head.updateJoint(0,lastJoint)
            result.head.joint(1) = null
          }
          last = (segment.leftEx, segment.left)
          lastJoint = segment(0)
        }
        if (segment(1)!=null){
          result ::= new Segment(f"${kv._1}$last${(segment.right, segment.rightEx)}")
          result.head.updateJoint(0,lastJoint)
          result.head.updateJoint(1,segment(1))
          last = (segment.right,segment.rightEx)
          lastJoint = null
        }
      }
      val tiglen = contigs(kv._1).length
      if (last._1!=tiglen) {
        result ::= new Segment(f"${kv._1}$last${(tiglen,tiglen)}")
        result.head.updateJoint(0, lastJoint)
        result.head.joint(1) = null
      }
      result
    }}
    println(f"[Info] ${segmentsByContig.size} segments from ${segmentsByContig.map(_.contig).toSet.size} contigs are involved in correction/scaffolding.")
    segmentsByContig
  }
  def mergeAndWrite(segments: Iterable[Segment]) = {
    val visit = mutable.Set[Segment]()
    def traverse(os:OrientedSegment, contig:Contig){
      if (visit.contains(os.segment)) return
      visit.add(os.segment)
      if (os.bestR!=null){
        val (link,idx) = os.bestR
        os.pruneR(link.overlap(1+idx)-1)
        link(idx^1,1).pruneL(link.overlap(3-idx*3))
        contig.append(os.getName,os.seq)
        traverse(link(idx^1,1), contig)
      } else contig.append(os.getName,os.seq)
    }
    for ( segment <- segments if segment.joint.contains(null) && !visit.contains(segment)) { // for chain
      val contig = new Contig()
      traverse(new OrientedSegment(segment,segment(0)==null), contig)
      outFasta.writeRecord(contig.name.toString, contig.buffer.toString)
    }
    for ( segment <- segments if !visit.contains(segment)) { // for circle
      val contig = new Contig()
      traverse(new OrientedSegment(segment,true), contig)
      outFasta.writeRecord(contig.name.toString, contig.buffer.toString)
    }
    segments.foreach(s=>contigs.remove(s.contig))
    contigs.foreach(kv=>outFasta.writeRecord(kv._1,kv._2))
  }

  def correctAndScaffoldFasta() = {
    loadFasta()
    val filteredSegments = linkFilter()
    mergeAndWrite(filteredSegments)
    outFasta.writer.close()
  }
  def put(a:String,b:String,s:Double) = {
    records ::= (a,b,s)
  }
  def loadEdgeFromFile(filename:String) = {
    for (line <- scala.io.Source.fromFile(filename).getLines() if line.startsWith("+") || line.startsWith("-")){
      val tokens = line.split(raw"\s+")
      records ::= (tokens(0),tokens(1),tokens(2).toDouble)
    }
  }
}
