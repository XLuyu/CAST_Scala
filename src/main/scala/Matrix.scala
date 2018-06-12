/**
  * Created by workshop on 14-Apr-18.
  */
class Genotype(var acgt:Array[Double]){
  val sum = acgt.sum
  val max = acgt.max
  def isUnreliable = sum < 5
  def proportioning() = {
    if (sum!=0) acgt = acgt.map(_/sum)
    this
  }
  def normalize() = {
    if (max!=0) acgt = acgt.map(_/max)
    this
  }
  def distance(other:Genotype) = {
    if (this.isUnreliable || other.isUnreliable)
      10.0
    else
      (this.acgt zip other.acgt map {case (x,y)=>Math.abs(x-y)}).sum / 2
  }
}

class GenotypeVector(val vector:Array[Genotype]){
  def toDistMatrix = {
      for (i <- vector) yield
        for (j <- vector) yield {
          val d = i.distance(j)
          if (d < 0.1) 0 else d
        }
  }
  def isReliable = {
    !vector.exists(_.isUnreliable)
  }
  def isHeterogeneousPrecheck = {
    var max = (0.0,0.0,0.0,0.0)
    var min = (1.0,1.0,1.0,1.0)
    for ( gt <- vector) {
      max = (Math.max(max._1,gt.acgt(0)),Math.max(max._2,gt.acgt(1)),Math.max(max._3,gt.acgt(2)),Math.max(max._4,gt.acgt(3)))
      min = (Math.min(min._1,gt.acgt(0)),Math.min(min._2,gt.acgt(1)),Math.min(min._3,gt.acgt(2)),Math.min(min._4,gt.acgt(3)))
    }
    (max._1-min._1)+(max._2-min._2)+(max._3-min._3)+(max._4-min._4)>=0.4
  }
  def isHeterogeneous:Boolean = {
    if (!isHeterogeneousPrecheck) return false
    !toDistMatrix.exists(x=>{ val sx = x.sorted; sx.last-sx(1)<0.2}) // to check if heterogeneity is not from mutation
  }
  def consistentWithDepthLimit(limit:Array[Double]) = {
     vector.indices.exists(i=>vector(i).sum>limit(i))
  }
}

class Matrix(var data:Array[Array[Double]]=null, size:Int=0){
  if (data==null) data = Array.ofDim[Double](size,size)
  def *=(factor:Double) = {
    for (i <- data.indices; j <- data(i).indices)
      data(i)(j) *= factor
    this
  }
  def +=(other:Matrix) = {
    for (i <- data.indices; j <- data(i).indices)
      data(i)(j) += other.data(i)(j)
    this
  }
  def copy = {
    val other = new Matrix(size=data.length)
    for ( i <- data.indices ; j <- data(i).indices)
      other.data(i)(j) = data(i)(j)
    other
  }
  def normalized = {
    val max = data.map(_.max)
    if (!max.contains(0.0)) {
      for ( i <- data.indices ; j <- data(i).indices)
        data(i)(j) /= max(i)
    }
    this
  }
  def printline(hint:String): Unit ={
    println(s"======$hint")
    for ( i <- data){
      i.foreach(j=>print(f"$j%.2f\t"))
      println()
    }
  }
}
