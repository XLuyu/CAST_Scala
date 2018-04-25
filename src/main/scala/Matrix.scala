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

class Matrix(size:Int, initData:Array[Array[Double]]=null){
  var data = if (initData!=null) initData else Array.ofDim[Double](size,size)
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
    val other = new Matrix(data.length)
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
