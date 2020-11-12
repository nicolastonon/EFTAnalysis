#ifndef common_TFModel_h
#define common_TFModel_h

#include <string>
#include <vector>

#include "CMSSW_TensorFlow.h" //LOCAL
//#include "PhysicsTools/TensorFlow/interface/TensorFlow.h" //CMSSW

class TFModel {
public:
    // Analysis routine
    TFModel(const std::string &model_name, const unsigned _n_inputs,
            const std::string &_input_name, const unsigned _n_outputs,
            const std::string &_output_name);
    ~TFModel();
    std::vector<float> evaluate(float inputs[]);
    std::vector<float> evaluate(std::vector<float>&); //Overload
    std::vector<float> evaluate(std::vector<float*>&); //Overload
    void evaluate_fast(tensorflow::Tensor&, std::vector<tensorflow::Tensor>&); //Overload

private:
    //-- CHANGED
    tensorflow::Session* session = nullptr;
    tensorflow::GraphDef* graphDef = nullptr;

    const unsigned n_inputs;
    const unsigned n_outputs;
    const std::string input_name;
    const std::string output_name;
};

#endif // ifndef common_Analyzer_h
