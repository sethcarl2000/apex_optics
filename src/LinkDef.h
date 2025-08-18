#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class NPoly+; 
#pragma link C++ struct NPoly::NPolyElem+; 
#pragma link C++ class RMatrix+; 
#pragma link C++ class NPolyArray+; 
#pragma link C++ class MultiLayerPerceptron+; 
#pragma link C++ struct MultiLayerPerceptron::WeightGradient_t+; 
#pragma link C++ struct MultiLayerPerceptron::HessianTensor_t+; 
#pragma link C++ class PolynomialCut+; 
#pragma link C++ struct PolynomialCut::Segment_t+; 
#pragma link C++ struct PolynomialCut::Vertex_t+; 
#pragma link C++ class PolynomialCut::InvalidVertexException+; 
#pragma link C++ class PolynomialCut::DBFileException+; 

#endif 