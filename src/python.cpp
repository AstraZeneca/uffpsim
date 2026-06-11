// @author Rajendra Kumar

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include "searchEngine.h"
#include "pythonDocString.h"

namespace py = pybind11;

PYBIND11_MODULE(uffpsimLib, m) {
    static py::exception<FileNotFoundException> fileNotFoundEx(m, "FileNotFoundError", PyExc_FileNotFoundError);
    py::register_exception_translator([](std::exception_ptr p) {
        try { if (p) std::rethrow_exception(p); }
        catch (const FileNotFoundException &e) { fileNotFoundEx(e.what()); }
    });

    m.doc() = R"pbdoc(
        uffpsimLib
        ---------
        .. currentmodule:: uffpsimLib
        .. autosummary::
           :toctree: _generate
    )pbdoc";

    m.def("getCompactFingerPrintArray", &utils::getCompactFingerPrintArray);
    m.def("doBitwiseAndGetOnesCount", &utils::bitwiseAndPopcount);

    py::class_<FingerprintStore>(m, "FingerprintStore")
        .def(py::init<std::string, int,std::string,int,std::string,std::string,float,std::string,bool>(), py::arg("filename"), py::arg("mol_id_max_chars")=15, py::arg("mode")="r", 
                                                       py::arg("fpSize")=1024, py::arg("fp_params")="", py::arg("info")="", py::arg("cluster_threshold")=0.15,
                                                       py::arg("cluster_mode")="memory", py::arg("cluster_parallel")=false)
        .def_property("fp_params_json", &FingerprintStore::getFingerprintParameters, nullptr, "input fingerprint parameters as JSON format.")
        .def_property("info", &FingerprintStore::getInfo, nullptr, "any additional information stored in the database file.")
        .def_property("mol_id_max_chars", &FingerprintStore::getMolIdMaxLength, nullptr, "maximum allowed length of molecule ID.")
        .def_property("fp_bits_size", &FingerprintStore::getFPBitsSize, nullptr, "size of fingerprint in bits.")
        .def_property("inner_clustering_threshold", &FingerprintStore::getInnerClusteringThreshold, nullptr, "inner clustering threshold.")
        .def_property("num_mols", &FingerprintStore::getNumberOfMolecules, nullptr, "total number of molecules in the database.")
        .def("set_info", &FingerprintStore::setInfo, FpStoreSetInfoDoc)
        .def("append_fingerprints", &FingerprintStore::AppendFingerprints, FpStoreAppendFingerprintsDoc)
        .def("perform_clustering_write", &FingerprintStore::performInnerClusteringAndWrite, FpStorePerformInnerClusteringAndWriteDoc, py::call_guard<py::gil_scoped_release>())
        .def("redo_clustering_write", &FingerprintStore::reDoInnerClusteringInMemory, FpStoreReDoInnerClusteringInMemory, py::call_guard<py::gil_scoped_release>())
        .def("magic_number_exists", &FingerprintStore::magicNumberExists, FpStoreMagicNumberExistsDoc)
        .def("set_magic_number", &FingerprintStore::setMagicNumber, FpStoreSetMagicNumberDoc)
        .def("get_smiles_for_id", &FingerprintStore::getSmilesFromID, FpStoreGetSmilesFromIdDoc)
        .def("build_mol_id_to_index_map", &FingerprintStore::buildMolIdToIndexMap, FpStorebuildMolIdToIndexMap)
        .def("build_mol_id_index_table", &FingerprintStore::buildMolIdIndexTable, FpStoreBuildMolIdIndexTableDoc)
        .def("update_fingerprints", &FingerprintStore::updateFingerprints, FpStoreUpdateFingerprintsDoc)
        .def("load_data_in_memory", &FingerprintStore::loadDataInMemory, FpStoreLoadDataInMemoryDoc)
        .def("close", &FingerprintStore::close, FpStoreCloseDoc);

    py::class_<FPSearchEngine>(m, "FPSearchEngineBase")
        .def(py::init<std::string, std::string>(), py::arg("filename"), py::arg("mode") = "memory", py::call_guard<py::gil_scoped_release>())
        .def_readonly("fp_store", &FPSearchEngine::_fpStore, "FingerprintStore object")
        .def_property("num_mols", &FPSearchEngine::getNumberOfMolecules, nullptr, "total number of molecules in the database.")
        .def("close", &FPSearchEngine::close, FpSearchEngineCloseDoc)
        .def("_search", &FPSearchEngine::search, py::call_guard<py::gil_scoped_release>(), FpSearchEngineSearchDoc)
        .def("_batch_search", &FPSearchEngine::batchSearch, py::call_guard<py::gil_scoped_release>(), FpSearchEngineBatchSearchDoc);
};
