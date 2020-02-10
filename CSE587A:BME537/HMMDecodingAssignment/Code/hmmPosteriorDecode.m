(* Wolfram Language package *)

(* posteriorDecode should return the state names for the sequence of most likely state.*)
posteriorDecode[observationSeq_, hmm_] :=
 Module[{posteriorMatrix,i},
 posteriorMatrix = posteriorProbabilities[observationSeq,hmm];
 
  (*Calculate most likely state and map it to the state representation in the hmm*)
 Table[hmm["states"][[Ordering[posteriorMatrix[[i]],-1][[1]]]],{i,Length[observationSeq]}]
 
 
 ]
posteriorProbabilities[observationSeq_, hmm_] := 
 Module[{forwardMatrix,backwardMatrix,posteriorMatrix,i},
 	
 forwardMatrix = buildForwardMatrix[observationSeq,hmm];
 backwardMatrix = buildBackwardMatrix[observationSeq,hmm];
 
 posteriorMatrix = forwardMatrix*backwardMatrix;
 Do[norm[posteriorMatrix[[i]]],{i,Length[observationSeq]}];
 
 posteriorMatrix
 
 
 ]
 buildForwardMatrix[observationSeq_, hmm_] := 
	Module[{numberOfObservations, numberOfStates, forwardMatrix, observationIndex,i,j,k,possibilities},
	  (* Put your code for bulding the forward matrix here. To give you an idea of what to expect,
         my code is 13 lines. But use as many lines as you need to make your code clear and readable. *)
	  numberOfObservations = Length[observationSeq];
	  numberOfStates = Length[hmm["states"]];
	  forwardMatrix = Table[0,{numberOfObservations},{numberOfStates}];
	  Do[
		
	forwardMatrix[[1,i]] = hmm["initialStateProbs"][[i]]*hmm["emissionMatrix"][[observationSeq[[1]],i]];
	,{i,numberOfStates}];
	 forwardMatrix[[1]] = norm[forwardMatrix[[1]]];
	 
	Do[
		
		Do[
			
			possibilities = Table[hmm["transitionMatrix"][[k,i]]*forwardMatrix[[j-1,k]],{k,numberOfStates}];
			
			observationIndex = observationSeq[[j]];
			forwardMatrix[[j,i]] = Total[hmm["emissionMatrix"][[observationIndex,i]]*possibilities];
			
			,{i,numberOfStates}];
			
		forwardMatrix[[j]] = norm[forwardMatrix[[j]]];
		
		,{j,Range[2,numberOfObservations]}]; 
	  
	  
	  
	 (* Return the Viterbi matrix. *)
	 forwardMatrix] 

buildBackwardMatrix[observationSeq_, hmm_] := 
	Module[{numberOfObservations, numberOfStates, backwardMatrix,i,j,k,possibilities, observationIndex},
	  (* Put your code for bulding the Viterbi matrix here. To give you an idea of what to expect,
         my code is 13 lines. But use as many lines as you need to make your code clear and readable. *)
 	  numberOfObservations = Length[observationSeq];
	  numberOfStates = Length[hmm["states"]];
	  backwardMatrix = Table[0,{numberOfObservations},{numberOfStates}];
	  Do[
		
	backwardMatrix[[numberOfObservations,i]] = 1./numberOfStates;
	,{i,numberOfStates}];
	 
	 
	Do[
		
		Do[
			observationIndex = observationSeq[[j+1]];
			possibilities = Table[hmm["transitionMatrix"][[i,k]]*hmm["emissionMatrix"][[observationIndex,k]]*backwardMatrix[[j+1,k]],{k,numberOfStates}];
			
			
			backwardMatrix[[j,i]] = Total[possibilities];
			
			,{i,numberOfStates}];
			
		backwardMatrix[[j]] = norm[backwardMatrix[[j]]];
		
		,{j,numberOfObservations-1,1,-1}]; 
	 (* Return the Viterbi matrix. *)
	 backwardMatrix] 
 	 
norm[viterbi_] := 
	Module[{},
		If[Total[viterbi]==0,viterbi/Total[viterbi],Normalize[viterbi,Total]]
	]