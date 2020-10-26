#include "TFModel.h"

TFModel::TFModel(const std::string &model_name, const unsigned _n_inputs,
    const std::string &_input_name, const unsigned _n_outputs,
    const std::string &_output_name):
    n_inputs(_n_inputs),
    n_outputs(_n_outputs),
    input_name(_input_name),
    output_name(_output_name)
{
    std::cout<<"Load tensorflow graph from "<<model_name<<std::endl;
    graphDef = tensorflow::loadGraphDef(model_name);

    //-- could also call sub-functions directly
    // graphDef = new tensorflow::GraphDef();
    // tensorflow::Status status = ReadBinaryProto(tensorflow::Env::Default(), model_name, graphDef);

    // if (!status.ok()) {
    //   std::cout<< "error while loading graph def: " << status.ToString()<<std::endl;
    // }

    // std::cout<<"Create tensorflow session"<<std::endl<<std::endl;

    // NICOLAS -- CHANGED -- GOT UNEXPLAINED ERROR : 'No session factory registered for the given session options'... => Get TF session directly
    // session = tensorflow::createSession((tensorflow::GraphDef*) graphDef); //1 thread by default

    //--> call sub-functions directly
    tensorflow::Status status;
    tensorflow::SessionOptions sessionOptions;
    status = NewSession(sessionOptions, &session);
    status = session->Create(*graphDef);
}

// std::vector<float> TFModel::evaluate(const double inputs[])
std::vector<float> TFModel::evaluate(float inputs[])
{
    tensorflow::Tensor input(tensorflow::DT_FLOAT, { 1, n_inputs });
    float* d = input.flat<float>().data();
    for(unsigned i=0; i < n_inputs; i++){
        *d++ = (float)inputs[i];
    }
    std::vector<tensorflow::Tensor> outputs;
    tensorflow::run((tensorflow::Session*)session, { { input_name, input } }, { output_name }, &outputs);

    std::vector<float> out;
    for(unsigned i=0; i < n_outputs; i++){
          out.push_back(outputs[0].matrix<float>()(0,i));
    }
    return out;
}

//Overload
std::vector<float> TFModel::evaluate(std::vector<float> v_inputs)
{
    tensorflow::Tensor input(tensorflow::DT_FLOAT, {1, n_inputs} );
    float* d = input.flat<float>().data();

    for(unsigned i=0; i < n_inputs; i++) {*d++ = (float) v_inputs[i];}

    std::vector<tensorflow::Tensor> outputs;
    tensorflow::run((tensorflow::Session*)session, { {input_name, input} }, {output_name}, &outputs);

    std::vector<float> out;
    for(unsigned i=0; i < n_outputs; i++) {out.push_back(outputs[0].matrix<float>()(0,i));}

    return out;
}

//Overload //Avoid un-necessary copies, pass all by reference (not much faster...)
void TFModel::evaluate_fast(tensorflow::Tensor& input, std::vector<tensorflow::Tensor>& outputs)
// void TFModel::evaluate_fast(tensorflow::Tensor input, std::vector<float>& out)
{
    // tensorflow::Tensor input(tensorflow::DT_FLOAT, {1, n_inputs} );
    // float* d = input.flat<float>().data();
    // for(unsigned i=0; i < n_inputs; i++) {*d++ = (float) v_inputs[i];}

    // std::vector<tensorflow::Tensor> outputs;
    tensorflow::run((tensorflow::Session*)session, { {input_name, input} }, {output_name}, &outputs);

    // std::vector<float> out;
    // for(unsigned i=0; i < n_outputs; i++) {out.push_back(outputs[0].matrix<float>()(0,i));}

    return;
}

TFModel::~TFModel()
{
    // std::cout << "~TFModel()"<<std::endl;

    if(session != nullptr)
    {
        // std::cout << "Close tensorflow session"<<std::endl;
        tensorflow::Session* s = (tensorflow::Session*) session;
        tensorflow::closeSession(s);
    }

    if(graphDef != nullptr)
    {
        // std::cout << "Delete tensorflow graph"<<std::endl;
        delete (tensorflow::GraphDef*) graphDef;
    }
}
