/*
+--------------------------------------------------------+
¦ GenA 1.0 - General Purpose Genetic Algorithm           ¦
+--------------------------------------------------------¦
¦ Copyright © 2013 Roman Chirikov (http://romashka.ch)   ¦
+--------------------------------------------------------¦
¦ Licensed under the MIT license.                        ¦
+--------------------------------------------------------+
*/

function GenA(fitnessFunction, variables, param)
{
	if(!fitnessFunction || typeof fitnessFunction != "function")
		throw "Error: fitness function is not defined or is not a function.";
	if(fitnessFunction.length < 1)
		throw "Error: fitness function does not accept arguments.";
	if(!variables)
		throw "Error: variables are not defined.";
	
	this.userFunction = fitnessFunction;
	
	if(!param) var param = {};

	var defaults = {
		expected: 1,
		precision: 0.01,
		iterations: 0,
		populationSize: 10,
		crossoverProb: 0.7,
		mutationProb: 0.2,
		bitMutationProb : 0.01,
		elitism: true,
		timeout: 10,
		sleep: 0,
		raw: false,
		chromLength: 20,
		verbose: false,
		onstart: function(){return false;},
		onstop: function(){return false;},
		oniteration: function(){return false;}
	};

	// set parameters or their defaults
	for(var k in defaults)
	{
		if(typeof param[k] == "undefined")
			this[k] = defaults[k];
		else
			this[k] = param[k];
	}
	
	if(!this.raw)
	{
		this.chromLength = 0;
		
		for(var i = 0; i < variables.length; i++)
		{
			var ticks = (variables[i].to - variables[i].from) / variables[i].step;
			var bits = Math.ceil( Math.log(ticks) / Math.log(2) );
			if(Math.pow(2, bits) == ticks) bits++;

			this.chromLength += bits;
			variables[i].numbits = bits;
			var maxbin = Math.pow(2, bits) - 1;
			var coef = maxbin / ticks;
			variables[i].coef = coef;
		}
		
		if(this.chromLength < 1)
			throw "Error: variables are not defined properly.";
		
		this.variables = variables;
	}
	else
	{
		if(this.chromLength < 1)
			throw "Error: chromosome length is not defined properly.";
		
		this.variables = [];
	}
	
	this.population = [];
	this.bestfitnessever = 0;

	this.run = function()
	{
		this.initPopulation();
		if(this.verbose)
		{
			console.log("> Initial population:");
			this.showpop();
		}

		this.forceStop = false;
		this.counter = 0;

		if(this.verbose) console.log("> Start main loop");
		this.onstart();
		
		this.timestart = (new Date()).getTime();
		this.mainLoop();
	}
	
	this.mainLoop = function()
	{
		if(this.iterations != 0 && this.counter >= this.iterations) {this.mainStop(); return true;} // required iterations done
		if(this.iterations == 0 && 
			!(this.bestfitnessever < this.expected - this.precision || 
			  this.bestfitnessever > this.expected + this.precision)
		) {this.mainStop(); return true;} // required accuracy achieved
		if(this.forceStop) {this.mainStop(); return true;} // manual stop
		if(this.iterations == 0 && this.timeout > 0 && (new Date()).getTime() - this.timestart > this.timeout * 1000) {this.mainStop(); return true;} // timeout
		
		this.newGeneration();
		if(!this.raw)
			this.oniteration({iteration: this.counter, variables: this.getVars(this.population[this.currentBestChromoIndex].bits), fitness: this.population[this.currentBestChromoIndex].fitness});
		else
			this.oniteration({iteration: this.counter, bits: this.population[this.currentBestChromoIndex].bits, fitness: this.population[this.currentBestChromoIndex].fitness});
		this.counter++;
		var that = this;
		setTimeout(function(){that.mainLoop();}, this.sleep*1000);
	}
	
	this.mainStop = function()
	{
		if(this.verbose)
		{
			console.log("> Finished. Final population:");
			this.showpop();
			console.log("> Iterations done: " + this.counter);
			console.log("> Best fitness: " + this.bestfitnessever);
			console.log("> Solution:");
			if(this.raw)
			{
				this.showChrom(this.bestchromoever);
			}
			else
			{
				var vars = this.getVars(this.bestchromoever.bits);
				for(var i = 0; i < this.variables.length; i++)
				{
					console.log("var "+(i+1)+" = "+vars[i]);
				}
			}
		}
		
		this.onstop({iteration: this.counter, variables: this.getVars(this.bestchromoever.bits), fitness: this.bestfitnessever});
	}

	this.stop = function()
	{
		this.forceStop = true;
	}

	this.fitFunction = function(bits)
	{
		if(this.raw)
		{
			return this.userFunction.apply(this, [bits]);
		}
		else
		{
			return this.userFunction.apply(this, this.getVars(bits));
		}
	}
	
	this.getVars = function(bits)
	{
		var vars = [];
		var shift = 0;

		for(var i = 0; i < this.variables.length; i++)
		{
			vars[i] = this.getValue(bits.slice(shift, shift + this.variables[i].numbits), this.variables[i]);
			shift += this.variables[i].numbits;
		}
		
		return vars;
	}

	this.getValue = function(bits, variable)
	{
		var raw = parseInt(bits.join(''), 2);
		raw = raw / variable.coef;
		raw = raw * variable.step;
		raw = raw + variable.from;
		raw  = Math.round(raw / variable.step) * variable.step;
		return raw;
	}

	this.showpop = function(population)
	{
		if(typeof population == "undefined") var population = this.population;
		for(var c = 0; c < population.length; c++) // for each chromosome
		{
			this.showChrom(population[c], c);
		}
		console.log("");
	}

	this.showChrom = function(chromosome, n)
	{
		if(typeof n == "undefined") var n = "";
		console.log("chrom "+n+": "+chromosome.bits.join("")+"; fitness: "+chromosome.fitness);
		//(hash "+parseInt(chromosome.bits.join(""), 2)+")
	}

	this.initPopulation = function()
	{
		for(var c = 0; c < this.populationSize; c++) // for each chromosome
		{
			this.population[c] = {bits: [], fitness: 0};
			for(var b = 0; b < this.chromLength; b++) // for each bit
			{
				this.population[c].bits[b] = this.binRandom();
			}
			this.population[c].fitness = this.fitFunction(this.population[c].bits);
		}
		
		this.bestfitnessever = this.bestFitness(this.population);
		this.currentBestChromoIndex = this.bestChromoIndex(this.population);
		this.bestchromoever = {bits: this.population[this.currentBestChromoIndex].bits.slice(), fitness: this.population[this.currentBestChromoIndex].fitness};
	}

	this.newGeneration = function()
	{
		var rand;
		var newGen = [];
		
		if(this.verbose) console.log("> Main loop iteration "+this.counter);

		var totalFitness = 0;
		for(var c = 0; c < this.populationSize; c++) // for each chromosome
		{
			totalFitness += this.population[c].fitness;
		}
		if(totalFitness == 0) totalFitness = 1;

		var coefficient = 1000 / totalFitness;
		
		if(this.verbose) console.log("> Roulette-wheel pair selection:");
		for(p = 0; p < this.populationSize/2; p++) // for each pair
		{
			//make pairs
			if(this.verbose) console.log("> Pair #"+p+":");
			for(var m = 0; m < 2; m++) // for each member
			{
				rand = this.random1000();
				var sum = 0;
				for(var c = 0; c < this.populationSize; c++)
				{
					sum += coefficient * this.population[c].fitness;
					if(rand <= sum)
					{
						newGen[p * 2 + m] = {bits: this.population[c].bits.slice(), fitness: this.population[c].fitness};
						break;
					}
				}
			}
			if(this.verbose)
			{
				this.showChrom(newGen[p * 2], p*2);
				this.showChrom(newGen[p * 2 + 1], p*2+1);
				console.log("");
			}
		}
		
		for(p = 0; p < this.populationSize/2; p++) // for each pair
		{
			c1 = 2 * p;
			c2 = 2 * p + 1;
			// crossover?
			if(this.random() < this.crossoverProb)
			{
				point = Math.floor(this.random() * 1000000) % (this.chromLength - 2);
				if(this.verbose) console.log("> Crossover pair "+p+" at point: "+point);
				for(b = point + 1; b < this.chromLength; b++) // exchange bits
				{
					tempbit = newGen[c1].bits[b];
					newGen[c1].bits[b] = newGen[c2].bits[b];
					newGen[c2].bits[b] = tempbit;
				}
				if(this.verbose)
				{
					console.log("> Result:");
					this.showChrom(newGen[c1], c1);
					this.showChrom(newGen[c2], c2);
				}
			}
		}
		
		for(var c = 0; c < this.populationSize; c++) // for each chromosome
		{
			// mutation?
			if(this.random() < this.mutationProb) // yes
			{
				if(this.verbose) console.log("> Mutation of chromosome "+c+":");
				for(b = 0; b < this.chromLength; b++) // for each bit
				{
					if(this.random() < this.bitMutationProb)
					{
						newGen[c].bits[b] = 1 - newGen[c].bits[b];
					}
				}
				if(this.verbose) console.log("chrom 0: "+newGen[c].bits.join(""));
			}

			newGen[c].fitness = this.fitFunction(newGen[c].bits);
		}
		
		if(this.elitism)
		{
			// elitism
			if(this.bestFitness(this.population) > this.bestFitness(newGen))
			{
				if(this.verbose) console.log("> Elitism occured");
				var bestchromoindex = this.bestChromoIndex(this.population);
				newGen[this.worstChromoIndex(newGen)] = {bits: this.population[bestchromoindex].bits.slice(), fitness: this.population[bestchromoindex].fitness};
			}
		}
		
		// replace generation
		for(var c = 0; c < this.populationSize; c++) // for each chromosome
		{
			this.population[c] = {bits: newGen[c].bits.slice(), fitness: newGen[c].fitness};
		}
		
		this.currentBestChromoIndex = this.bestChromoIndex(this.population);
		
		if(this.bestFitness(this.population) > this.bestfitnessever)
		{
			this.bestfitnessever = this.bestFitness(this.population);
			this.bestchromoever = {bits: this.population[this.currentBestChromoIndex].bits.slice(), fitness: this.population[this.currentBestChromoIndex].fitness};
		}
		
		if(this.verbose)
		{
			console.log("> New generation:");
			this.showpop();
		}
	}

	this.bestFitness = function(population)
	{
		var bestFitness = -9999;
		for(var c = 0; c < this.populationSize; c++)
		{
			if(population[c].fitness > bestFitness)
			{
				bestFitness = population[c].fitness;
			}
		}
		return bestFitness;
	}

	this.bestChromoIndex = function(population)
	{
		var bestIndex = 0;
		var bestFitness = -9999;
		for(var c = 0; c < this.populationSize; c++)
		{
			if(population[c].fitness > bestFitness)
			{
				bestFitness = population[c].fitness;
				bestIndex = c;
			}
		}
		return bestIndex;
	}

	this.worstChromoIndex = function(population)
	{
		var bestIndex = 0;
		var bestFitness = 999999;
		for(var c = 0; c < this.populationSize; c++)
		{
			if(population[c].fitness < bestFitness)
			{
				bestFitness = population[c].fitness;
				bestIndex = c;
			}
		}
		return bestIndex;
	}

	this.random = Math.random;

	this.random1000 = function()
	{
		return Math.floor(Math.random() * 1000);
	}

	this.binRandom = function()
	{
		return Math.floor(Math.random() * 1000) % 2;
	}

	return this;
}