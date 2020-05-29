#For later prediction using C++ API, it is necessary to freeze the graph and make it static (see online docs)

import tensorflow
from Utils.ColoredPrintout import colors
from tensorflow.keras.models import load_model

# //--------------------------------------------
# //--------------------------------------------

def freeze_session(session, keep_var_names=None, output_names=None, clear_devices=True):
    """
    Freezes the state of a session into a pruned computation graph.

    Creates a new computation graph where variable nodes are replaced by
    constants taking their current value in the session. The new graph will be
    pruned so subgraphs that are not necessary to compute the requested
    outputs are removed.
    @param session The TensorFlow session to be frozen.
    @param keep_var_names A list of variable names that should not be frozen,
                          or None to freeze all the variables in the graph.
    @param output_names Names of the relevant graph outputs.
    @param clear_devices Remove the device directives from the graph for better portability.
    @return The frozen graph definition.
    """
    graph = session.graph
    with graph.as_default():
        freeze_var_names = list(set(v.op.name for v in tensorflow.compat.v1.global_variables()).difference(keep_var_names or []))
        output_names = output_names or []
        output_names += [v.op.name for v in tensorflow.compat.v1.global_variables()]
        graphdef_inf = tensorflow.compat.v1.graph_util.remove_training_nodes(graph.as_graph_def())
        if clear_devices:
            for node in graphdef_inf.node:
                node.device = ""
        frozen_graph = tensorflow.compat.v1.graph_util.convert_variables_to_constants(
            session, graphdef_inf, output_names, freeze_var_names)
        return frozen_graph

# //--------------------------------------------
# //--------------------------------------------

def FreezeSession_and_SaveModel(opts, sess, weightDir, h5modelName):

    print('\n'); print(colors.fg.lightblue, "--- Freeze graph...", colors.reset); print('\n')

    tensorflow.keras.backend.set_learning_phase(0) # This line must be executed before loading Keras model (else mismatch between training/eval layers, e.g. Dropout)
    # model = load_model(h5modelName) # model has to be re-loaded
    model = load_model(h5modelName, compile=False) # model has to be re-loaded #compile=False <-> does not need to define any custom loss, since not needed for testing

    # tensorflow.compat.v1.keras.backend.clear_session() #Closing the last session avoids that node names get a suffix appened when opening a new session #Does not work?
    # sess = tensorflow.compat.v1.keras.backend.get_session()
    # graph = sess.graph

    inputs_names = [input.op.name for input in model.inputs]
    outputs_names = [output.op.name for output in model.outputs]
    # print('\ninputs: ', model.inputs)
    print('\n')
    print(colors.fg.lightgrey, '--> inputs_names: ', inputs_names[0], colors.reset, '\n')
    # print('\noutputs: ', model.outputs)
    print(colors.fg.lightgrey, '--> outputs_names: ', outputs_names[0], colors.reset, '\n')
    # tf_node_list = [n.name for n in  tensorflow.compat.v1.get_default_graph().as_graph_def().node]; print('nodes list : ', tf_node_list)
    frozen_graph = freeze_session(sess, output_names=[output.op.name for output in model.outputs])
    tensorflow.io.write_graph(frozen_graph, weightDir, 'model.pbtxt', as_text=True)
    tensorflow.io.write_graph(frozen_graph, weightDir, 'model.pb', as_text=False)
    print('\n'); print(colors.fg.lightgrey, '===> Successfully froze graph :', colors.reset, weightDir+'model.pb', '\n')

    #-- Also append the names of the input/output nodes in the file "NN_info.txt" containing input features names, etc. (for later use in C++ code)
    text_file = open(weightDir + "NN_info.txt", "a+") #Append mode
    text_file.write(inputs_names[0]); text_file.write(' -1 -1 \n'); #use end values as flags to signal these lines
    text_file.write(outputs_names[0]); text_file.write(' -2 -2 \n');
    text_file.write(str(opts["nofOutputNodes"])); text_file.write(' -3 -3 \n');
    text_file.close()

    return

# //--------------------------------------------
# //--------------------------------------------
