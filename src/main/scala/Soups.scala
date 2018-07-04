import java.io.File
import org.tc33.jheatchart.HeatChart
import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer
import htsjdk._

class Soups(val bamFilenames: Array[String]) {
  val decay = 0.99995
  val coverage = new DepthDetector(bamFilenames).getCoverage
  val bamScanner = new BamFilesParallelScanner(bamFilenames)
  val stringer = new Stringer()
  print(coverage.map(_.toInt).mkString("\nDepth:","|","\n"))

  def support(a: Array[Array[Double]], b: Array[Array[Double]]) = {
    (for (i <- a.indices) yield {
      val q = a(i).indices.filter(_ != i).map(j => a(i)(j))
      val p = b(i).indices.filter(_ != i).map(j => b(i)(j))
      Math.abs(Util.PearsonCorrelationSimilarity(q, p))
    }).sum / a.length
  }
  def getOneChromGVFromScanner: ArrayBuffer[(String, (Matrix, Matrix))] = {
    // scan this contig to get all snp sites
    val contig = bamScanner.getChrom
    val contiglen = bamScanner.header.getSequence(contig).getSequenceLength
    var sites = new ArrayBuffer[(Int, GenotypeVector)](contiglen / 1000)
    do {
      val gv = new GenotypeVector(bamScanner.get())
      if (gv.isReliable && gv.isHeterogeneous) {
        sites.append((bamScanner.pos, gv))
      }
    } while (bamScanner.nextPosition())
    sites = sites.filter { case (_, gv) => gv.consistentWithDepth(coverage) }
    if (sites.length < 2) {
      println(f"\n[$contig] SNP: ${sites.length} sites (No enough SNP)")
      return new ArrayBuffer[(String, (Matrix, Matrix))]()
    } else {
      println(f"\n[$contig] SNP: ${sites.length} sites")
    }
    // train matrix
    var rightScanSum, leftScanSum, matrix = new Matrix(size = bamFilenames.length)
    val snapshot = new ArrayBuffer[Matrix](sites.length)
    val segment = new ArrayBuffer[(((Int, Int), (Int, Int)), (Matrix, Matrix))]()
    var lastpos = 0
    matrix.fill(Double.NaN)
    for ((pos, gv) <- sites) {
      matrix.updateBy(new Matrix(gv.toDistMatrix))
      if (!matrix.containsNaN) {
        val decayN = Math.pow(decay, pos - lastpos)
        leftScanSum *= decayN += (matrix * ((1 - decayN) / (1 - decay)))
        snapshot.append(leftScanSum.copy)
      } else
        snapshot.append(null)
      lastpos = pos
    }
    if (snapshot.last == null) {
      println(f"\n[$contig] SNP: some sample is missing across whole contig")
      for (i <- matrix.data.indices if matrix.data(i).forall(_.isNaN)) print(f"sample $i is missing\n")
      return new ArrayBuffer[(String, (Matrix, Matrix))]()
    }
    lastpos = contiglen + 1
    matrix.fill(Double.NaN)
    var lastSegEnd = sites.length - 1
    for (i <- sites.indices.reverse) {
      val (pos, gv) = sites(i)
      matrix.updateBy(new Matrix(gv.toDistMatrix))
      if (!matrix.containsNaN) {
        val decayN = Math.pow(decay, lastpos - pos)
        rightScanSum *= decayN += (matrix * ((1 - decayN) / (1 - decay)))
        if (0 < i && snapshot(i - 1) != null && support(snapshot(i - 1).data, rightScanSum.data) < 0.3) {
          val seg = (((sites(i - 1)._1, sites(i)._1), (
            if (lastSegEnd == sites.length - 1) contiglen else sites(lastSegEnd)._1,
            if (lastSegEnd == sites.length - 1) contiglen else sites(lastSegEnd + 1)._1)),
            (rightScanSum.copy, snapshot(lastSegEnd)))
          segment.append(seg)
          lastSegEnd = i - 1
        }
      }
      lastpos = pos
    }
    val seg = (((1, 1), (
      if (lastSegEnd == sites.length - 1) contiglen else sites(lastSegEnd)._1,
      if (lastSegEnd == sites.length - 1) contiglen else sites(lastSegEnd + 1)._1)),
      (rightScanSum.copy, snapshot(lastSegEnd)))
    segment.append(seg)
    for (i <- 1 until segment.length - 1) {
      val (j, k, l) = (segment(i - 1), segment(i), segment(i + 1))
      if (j._1._1 == k._1._2 && k._1._1 == l._1._2) {
        segment(i + 1) = ((l._1._1, k._1._2), (l._2._1, k._2._2))
        segment(i) = (k._1,null)
      }
    }
    val named_segment = segment.filter(_ != null).map(x => (f"$contig${x._1._1}${x._1._2}", x._2))
    for ((name, (headMatrix, tailMatrix)) <- named_segment) {
      new HeatChart(headMatrix.data).saveToFile(new File(s"png/${name}_L.png"))
      new HeatChart(tailMatrix.data).saveToFile(new File(s"png/${name}_R.png"))
    }
    named_segment
  }
  def getAllGVFromScanner = {
    val GVmap = mutable.Map[String, (Matrix, Matrix)]()
    while (bamScanner.hasNext) {
      val CGV = getOneChromGVFromScanner
      for ((contig, (headMatrix, tailMatrix)) <- CGV)
        GVmap(contig) = (headMatrix, tailMatrix)
    }
    GVmap
  }
  def pairwiseMutualBest(contigsGV: mutable.Map[String, Matrix]) = {
    val candidate = ArrayBuffer[(String, (String, Double))]()
    for ((k, v) <- contigsGV) yield {
      val contig = k.substring(1)
      val similarities = contigsGV.filter(_._1.substring(1) != contig).mapValues(v2 => support(v.data, v2.data))
      similarities.filter(_._2 > 0.6).foreach(x => candidate.append((k, x)))
    }
    candidate
  }
  def run(): Unit = {
    //    val coverages = bamScanner.getCoverage
    //    println(coverages.mkString(" "))
    val contigsTwoEndGV = getAllGVFromScanner
    val contigsOneEndGV = contigsTwoEndGV.flatMap { case (contig, (v1, v2)) => List(("+" + contig, v1), ("-" + contig, v2)) }
    val finalLink = pairwiseMutualBest(contigsOneEndGV)
    for ((k, v) <- finalLink if k < v._1) {
      println(f"$k\t${v._1}\t${v._2}")
      //      stringer.put(k,v._1,v._2)
    }
    //    stringer.getPath.foreach(println)
  }
}

object Soups {
  def main(args: Array[String]): Unit = new Soups(args).run()
}
