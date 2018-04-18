import scala.collection.mutable

class Stringer {
  val edge = mutable.Map[String,String]()
  val degree = mutable.Map[String,Int]().withDefaultValue(0)
  val weight = mutable.Map[String,Double]()
  def put(a:String,b:String,s:Double) = {
    edge(a) = b
    edge(b) = a
    weight(a) = s
    weight(b) = s
    degree(a.substring(1)) += 1
    degree(b.substring(1)) += 1
  }
  def getPath:mutable.Iterable[String] = {
    val visit = mutable.Set[String]()
    def dfs(endpoint:String):String = {
      visit.add(endpoint.substring(1))
      if (!edge.contains(endpoint)) return ""
      val target = edge(endpoint)
      val tcontig = target.substring(1)
      if (target(0)=='+')
        f"\t--${weight(endpoint)}%.2f--> "+tcontig+dfs("-"+tcontig)
      else
        f"\t--${weight(endpoint)}%.2f--> ~"+tcontig+dfs("+"+tcontig)
    }
    for ( (contig,d) <- degree if d==1 && !visit.contains(contig)) yield
      if (edge.contains("+"+contig))
        "~"+contig+dfs("+"+contig)
      else
        contig+dfs("-"+contig)
  }
}
