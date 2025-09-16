#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class NPoly+; 
#pragma link C++ struct SieveHole+; 
#pragma link C++ struct NPoly::NPolyElem+; 
#pragma link C++ struct ApexOptics::Trajectory_t+; 
#pragma link C++ struct ApexOptics::OpticsTarget_t+; 
#pragma link C++ class RMatrix+; 
#pragma link C++ class NPolyArray+; 
#pragma link C++ class NPolyArrayChain+; 
#pragma link C++ class MultiLayerPerceptron+; 
#pragma link C++ struct MultiLayerPerceptron::WeightGradient_t+; 
#pragma link C++ struct MultiLayerPerceptron::HessianTensor_t+; 

#pragma link C++ class PolynomialCut+; 
#pragma link C++ struct PolynomialCut::Segment_t+; 
#pragma link C++ struct PolynomialCut::Vertex_t+; 
#pragma link C++ class PolynomialCut::InvalidVertexException+; 
#pragma link C++ class PolynomialCut::DBFileException+; 
#pragma link C++ class PolynomialCutApp+; 

#pragma link C++ class PickSieveHoleApp+; 
#pragma link C++ class SaveOutputFrame+; 
#pragma link C++ class EvaluateCutFrame+; 
#pragma link C++ struct FPcoordPolynomial+; 
#pragma link C++ struct SieveHoleData+; 

#endif 