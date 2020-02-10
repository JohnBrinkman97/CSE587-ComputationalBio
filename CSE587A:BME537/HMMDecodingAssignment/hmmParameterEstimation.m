(* Wolfram Language package *)
(* Use this file to write code for supervised and unsupervised HMM parameter estimation.*)

(* supervisedHMMParameterEstimation should carry out maximum likelihood parameter estimation
   based on a set of observations and state labels, such as those provided in Test/mixed2.fa
   and Test/mixed2key.fa. It should output an HMM object that can pass checkHMMValidity. All
   the information needed to create an HMM object, including the number and names of the states
   and observation letters, should be taken from the two provided arguments. 
  *) 
supervisedHMMParameterEstimation[observations_, stateSequence_]:=
Module[{transitionMatrix,observationCounts,observationPercentages,stateCounts,stateProbs,states,HMMobject,emissionMatrix,emissionMatrixCounts,numberOfTransitions},
	
observationCounts = {0,0,0,0};
Map[observationCounts[[#]]++&,observations];
observationPercentages = observationCounts/Length[observations];
stateCounts = ConstantArray[0,CountDistinct[stateSequence]];
states = {};
Map[If[!MemberQ[states,#],AppendTo[states,#]]&,stateSequence];

Map[(stateCounts[[Position[states,#][[1]]]]++)&,stateSequence];

emissionMatrixCounts = ConstantArray[0,{Length[observationCounts],Length[states]}];
Do[emissionMatrixCounts[[observations[[i]],Position[states,stateSequence[[i]]][[1]]]]++,{i,Length[observations]}];
emissionMatrix = emissionMatrixCounts / Length[observations];
stateProbs = N[stateCounts/Length[stateSequence]];
numberOfTransitions = 0;
Map[If[stateSequence[[#-1]]!=stateSequence[[#]],numberOfTransitions++]&,Range[2,Length[stateSequence]]];
transitionMatrix = {};
	AppendTo[transitionMatrix,{(stateCounts[[1]]/numberOfTransitions/Length[states]-1)/(stateCounts[[1]]/numberOfTransitions/Length[states]),(stateCounts[[1]]/numberOfTransitions/Length[states])/(stateCounts[[1]]/numberOfTransitions/Length[states])}];
	AppendTo[transitionMatrix,{(stateCounts[[2]]/numberOfTransitions/Length[states])/(stateCounts[[2]]/numberOfTransitions/Length[states]),(stateCounts[[2]]/numberOfTransitions/Length[states]-1)/(stateCounts[[2]]/numberOfTransitions/Length[states])}];

HMMobject["states"] = states;
HMMobject["initialStateProbs"] = stateProbs;
HMMobject["transitionMatrix"] = transitionMatrix;
HMMobject["alphabet"] = {"A","C","G","T"};
HMMobject["emissionMatrix"] = emissionMatrix;
HMMobject

]
(* unsupervisedHMMParameterEstimation should carry out EM/ForwardBackward parameter estimation
   based on a set of observations. No state labels are provided. It should output an HMM object 
   that can pass checkHMMValidity. The number and names of the and observation letters, should 
   be taken from the first argument. The second argument is a list of strings that will serve
   as state names. It's length determines the number of states. 

unsupervisedHMMParameterEstimation[observations_, stateNames_]:=

*)