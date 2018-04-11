import java.io.File
import java.util

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
    val v2 = v.map(x=>if (x>=sum*0.9) 1 else 0 )
    if (v2.sum!=1 || sum<5) null else v2.map(_.toDouble)
//    if (sum<5) null else v.map(_/sum)
  }
  def distance(a:Array[Double],b:Array[Double]) = if (a==null || b==null) 10.0 else (a zip b map {case (x,y)=>(x-y)*(x-y)}).sum / 2

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
    if (GVSpanTable.nonEmpty && GVSpanTable.last.right-GVSpanTable.last.left<bamScanner.maxReadSpan) GVSpanTable.trimEnd(1)
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
  def getAllGVFromScanner ={
    val GVmap = mutable.Map[String,(Array[Int],Array[Int])]()
    while (bamScanner.hasNext) {
      val CGV = getOneChromGVFromScanner
      if (CGV._2!=null) GVmap(CGV._1) = CGV._2
    }
    GVmap
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
  def run(): Unit ={
//    val coverages = bamScanner.getCoverage
//    println(coverages.mkString(" "))
    val contigsTwoEndGV = getAllGVFromScanner
    val contigsOneEndGV = contigsTwoEndGV.flatMap {case (contig,(v1,v2)) => List(("+"+contig,v1),("-"+contig,v2))}
    val finalLink = pairwiseMutualBest(contigsOneEndGV)
    for ( (k,v) <- finalLink if k<v) {
      println(k,v)
      stringer.put(k,v)
    }
    stringer.getPath().foreach(println)
  }
}

object Soups {
  def main(args:Array[String]): Unit = new Soups(args).run()
}
