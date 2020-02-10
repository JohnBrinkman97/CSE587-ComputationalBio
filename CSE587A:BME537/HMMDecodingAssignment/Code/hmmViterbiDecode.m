(* decode
INPUT 
- observationSeq, in the form output by readFasta in tools.m
- hmm, in form output by readHMM in tools.m.
See tools.m for details.
     
OUTPUT
- stateSeq = a list containing the most likely state sequence, e.g., {h,h,m,m,m}.  
The length of this sequence equal to the length to the input observationSeq.  

IMPLEMENTATION
I recommend the following:
1) Implement the Viterbi table as a matrix that is transposed relative to the way
   it is shown in class and course notes. In other words, the rows correspond to
   observations and the columns to states.
2) After computing the Viterbi probabilities for each observation, normalize
   them by dividing by their sum. The matrix entries are no longer the Viterbi 
   probabilities, but the entries for each observation remain proportional
   to the Viterbi probabilities. The only thing you will use this matrix for is picking
   the state with the greatest entry for each observation, so scaling them all by a
   constant factor won't change the result. The benefit of doing this is that you
   avoid numerical underflow.
*)

viterbiDecode[observationSeq_, hmm_] :=
  Module[{numericStateSeq},
  (* Put your code here. *)
  (* Return the sequence of state names corresponding to the Viterbi decode path. *)
  
	Map[hmm["states"][[#]]&,traceback[buildMatrix[observationSeq,hmm],hmm]]
 ]

buildMatrix[observationSeq_, hmm_] := 
	Module[{numberOfObservations, numberOfStates, viterbiMatrix, observationIndex,i,j,k,possibilities},
	  (* Put your code for bulding the Viterbi matrix here. To give you an idea of what to expect,
         my code is 13 lines. But use as many lines as you need to make your code clear and readable. *)
	
	numberOfObservations = Length[observationSeq];
	numberOfStates = Length[hmm["states"]];
	viterbiMatrix = Table[0.,{numberOfObservations},{numberOfStates}];
	Do[
		
	viterbiMatrix[[1,i]] = hmm["initialStateProbs"][[i]]*hmm["emissionMatrix"][[observationSeq[[1]],i]];
	,{i,numberOfStates}];
	
	viterbiMatrix[[1]] = norm[viterbiMatrix[[1]]];
	Do[
		
		Do[
			
			possibilities = Table[hmm["transitionMatrix"][[k,i]]*viterbiMatrix[[j-1,k]],{k,numberOfStates}];
			
			observationIndex = observationSeq[[j]];
			viterbiMatrix[[j,i]] = hmm["emissionMatrix"][[observationIndex,i]]*Max[possibilities];
			
			,{i,numberOfStates}];
			
		viterbiMatrix[[j]] = norm[viterbiMatrix[[j]]];
		
		,{j,Range[2,numberOfObservations]}];
	 (* Return the Viterbi matrix. *)
	 viterbiMatrix] 
	 
	 
norm[viterbi_] := 
	Module[{},
		If[Total[viterbi]==0,viterbi/Total[viterbi],Normalize[viterbi,Total]]
	]


traceback[viterbiMatrix_, hmm_] :=
  (* Put your code for tracing back through the Viterbi matrix here. To give you an idea of what to expect,
     my code is 10 lines. But use as many lines as you need to make your code clear and readable. *)
	Module[{stateSeq,seqLength,numberOfStates,i,j,possibilities},
		seqLength = Length[viterbiMatrix];
		stateSeq = Table[0,seqLength];
	 
		numberOfStates = Length[hmm["states"]];
		stateSeq[[seqLength]] = Ordering[viterbiMatrix[[seqLength]],-1][[1]];
		
		Do[
		
			possibilities = Table[viterbiMatrix[[j,i]]*hmm["transitionMatrix"][[i,stateSeq[[j+1]]]],{i,numberOfStates}];
			
			stateSeq[[j]] = Min[Ordering[possibilities,-1]];
			
		,{j,seqLength-1,1,-1}];
		
	
	 (* Return the list of states. *)
	 stateSeq
	]	   
 
	