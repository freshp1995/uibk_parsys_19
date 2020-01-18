//pi program
use Random;

//var sample: int = 1000000000;
var sample: int = 100000;
var count: int = 0;

var seed = 13;

var randStream = new owned NPBRandomStream(real, seed, parSafe=false);

var sampleArray = {1..sample};


//create x y list with random numbers
count = + reduce [(x,y) in zip(randStream.iterate(sampleArray), randStream.iterate(sampleArray))] (x**2 + y**2) <= 1.0;


	

writeln(count * 4.0 / sample);
