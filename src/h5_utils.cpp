/**
 * @brief These are HDF5 utility functions used in this package.
 * @author Rajendra Kumar
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <H5Cpp.h>

#include "lib.h"

namespace h5 = H5;

std::string utils::randomString(size_t length)
{
    auto randchar = []() -> char
    {
        const char charset[] =
            "0123456789"
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            "abcdefghijklmnopqrstuvwxyz";
        const size_t max_index = (sizeof(charset) - 1);
        return charset[rand() % max_index];
    };
    std::string str(length, 0);
    std::generate_n(str.begin(), length, randchar);
    return str;
}

h5::Group *utils::createOrOpenGroup(h5::Group *parentGroup, const std::string &name)
{
    h5::Group *group;
    if (!parentGroup->nameExists(name))
    {
        group = new h5::Group(parentGroup->createGroup(name));
    }
    else
    {
        group = new h5::Group(parentGroup->openGroup(name));
    }
    return group;
}

h5::Group *utils::removeAndOpenGroup(h5::Group *parentGroup, const std::string &name)
{
    if (parentGroup->nameExists(name))
    {
        try
        {
            parentGroup->unlink(name);
        }
        catch (const H5::GroupIException &unlink_exception)
        {
            std::cerr << "Failed to unlink " << name << ". moving it to new random location." << std::endl;
            parentGroup->moveLink(name, name + "_toBeRemove_" + utils::randomString(8));
        }
    }
    return new h5::Group(parentGroup->createGroup(name));
}

template <typename T>
bool utils::set_scaler_attribute(h5::Group *group, const std::string &name, h5::PredType dtype, const T value)
{
    if (group->attrExists(name))
    {
        group->removeAttr(name);
    }
    h5::DataSpace attrSpace(H5S_SCALAR);
    h5::Attribute attr = group->createAttribute(name, dtype, attrSpace);
    attr.write(dtype, &value);
    attr.close();
    attrSpace.close();
    return true;
}
template bool utils::set_scaler_attribute<int>(h5::Group *group, const std::string &name, h5::PredType dtype, const int value);
template bool utils::set_scaler_attribute<float>(h5::Group *group, const std::string &name, h5::PredType dtype, const float value);
template bool utils::set_scaler_attribute<uint64_t>(h5::Group *group, const std::string &name, h5::PredType dtype, const uint64_t value);

template <typename T>
T utils::get_scaler_attribute(h5::Group *group, const std::string &name)
{
    if (!group->attrExists(name))
    {
        return 0;
    }

    h5::Attribute *attr = new h5::Attribute(group->openAttribute(name));
    T value;
    attr->read(attr->getDataType(), &value);
    attr->close();
    delete attr;
    return value;
}
template int utils::get_scaler_attribute<int>(h5::Group *group, const std::string &name);
template float utils::get_scaler_attribute<float>(h5::Group *group, const std::string &name);
template uint64_t utils::get_scaler_attribute<uint64_t>(h5::Group *group, const std::string &name);

bool utils::set_string_attribute(h5::Group *group, const std::string name, const std::string value)
{
    if (group->attrExists(name))
    {
        group->removeAttr(name);
    }
    h5::StrType attrType(h5::PredType::C_S1, H5T_VARIABLE);
    h5::DataSpace attrSpace(H5S_SCALAR);
    h5::Attribute magicFingerPrintWordAttrs = group->createAttribute(name, attrType, attrSpace);
    magicFingerPrintWordAttrs.write(attrType, value);

    magicFingerPrintWordAttrs.close();
    attrType.close();
    attrSpace.close();
    return true;
}

std::string utils::get_string_attribute(h5::Group *group, const std::string name)
{
    if (group->attrExists(name))
    {
        h5::Attribute *attr = new h5::Attribute(group->openAttribute(name));
        h5::StrType attrType = attr->getStrType();
        std::string result;
        attr->read(attrType, result);

        attr->close();
        delete attr;
        attrType.close();
        return result;
    }
    return "not-found";
}

void utils::addDataSetToH5(h5::Group *group, const std::string setName, h5::PredType dtype, const void *data, std::vector<size_t> dims)
{
    if (group->nameExists(setName))
    {
        try
        {
            group->unlink(setName);
        }
        catch (const H5::GroupIException &unlink_exception)
        {
            std::cerr << "Failed to unlink " << setName << ". moving it to new random location." << std::endl;
            group->moveLink(setName, setName + "_toBeRemove_-" + utils::randomString(8));
        }
    }

    hsize_t datasetDims[dims.size()];
    for (size_t i = 0; i < dims.size(); i++)
    {
        datasetDims[i] = dims[i];
    }
    h5::DataSpace *dataSpace = new h5::DataSpace(1, datasetDims);
    h5::DataSet *dataset = new h5::DataSet(group->createDataSet(setName, dtype, *dataSpace));
    dataset->write(data, dtype);

    dataSpace->close();
    dataset->close();
    delete dataset;
    delete dataSpace;
}

template <typename T>
bool utils::readDataSetFromH5(h5::Group *group, const std::string setName, h5::PredType dtype, std::vector<T> &output)
{
    if (!group->nameExists(setName)) {
        return false;
    }

    h5::DataSet *dataset = new h5::DataSet(group->openDataSet(setName));
    h5::DataSpace *dataSpace = new h5::DataSpace(dataset->getSpace());
    hsize_t dims_out[1];
    dataSpace->getSimpleExtentDims(dims_out, NULL);

    // std::cout << "dims_out[0] = " << dims_out[0] << std::endl;
    output.resize(dims_out[0]);
    dataset->read(output.data(), dataset->getDataType());

    dataSpace->close();
    dataset->close();
    delete dataSpace;
    delete dataset;
    
    return true;
}
template bool utils::readDataSetFromH5<int>(h5::Group*, const std::string, h5::PredType, std::vector<int>&);
template bool utils::readDataSetFromH5<uint64_t>(h5::Group*, const std::string, h5::PredType, std::vector<uint64_t>&);

void utils::appendFingerPrintDatasetToH5(h5::Group *group, const std::string setName, uint64_t *cfp, hsize_t size, hsize_t CFPSize)
{
    hsize_t offset[1] = {0};
    h5::PredType dtype = h5::PredType::NATIVE_UINT64;
    hsize_t maxdims[1]    = {H5S_UNLIMITED};
    if (!group->nameExists(setName)) {
        hsize_t chunk_dims[1] = {(hsize_t) 1000*CFPSize};
        

        // Modify dataset creation property to enable chunking
        h5::DSetCreatPropList prop;
        prop.setChunk(1, chunk_dims);

        hsize_t datasetDims[1] = {size};
        h5::DataSpace *dataSpace = new h5::DataSpace(1, datasetDims, maxdims);
        h5::DataSet *dataset = new h5::DataSet(group->createDataSet(setName, dtype, *dataSpace, prop));
        dataset->write(cfp, dtype);

        dataset->close();
        dataSpace->close();
        delete dataSpace;
        delete dataset;
    } else {
        //open dataset
        h5::DataSet *dataset = new h5::DataSet(group->openDataSet(setName));

        // get existing size from file
        h5::DataSpace *fileSpace = new h5::DataSpace(dataset->getSpace());
        fileSpace->getSimpleExtentDims(offset, NULL);
        delete fileSpace;

        // extend size of dataset
        hsize_t totalSize[1] = {offset[0] + size};
        dataset->extend(totalSize);
        
        // select hyperslab space
        hsize_t datasetDims[1] = {size};
        fileSpace = new h5::DataSpace(dataset->getSpace());
        fileSpace->selectHyperslab(H5S_SELECT_SET, datasetDims, offset);
        
        // write extended data
        h5::DataSpace *memSpace = new h5::DataSpace(1, datasetDims, NULL);
        dataset->write(cfp, dtype, *memSpace, *fileSpace);

        fileSpace->close();
        memSpace->close();
        dataset->close();
        delete fileSpace;
        delete memSpace;
        delete dataset;
    }
}

void utils::createEmptyFpDataset(h5::Group *group, const std::string setName, hsize_t size, hsize_t CFPSize)
{
    h5::PredType dtype = h5::PredType::NATIVE_UINT64;
    hsize_t maxdims[1]    = {H5S_UNLIMITED};
    if (!group->nameExists(setName)) {
        hsize_t chunk_dims[1] = {(hsize_t) 1000*CFPSize};
        

        // Modify dataset creation property to enable chunking
        h5::DSetCreatPropList prop;
        prop.setChunk(1, chunk_dims);
        prop.setFillValue(dtype, {0});

        hsize_t datasetDims[1] = {size};
        h5::DataSpace *dataSpace = new h5::DataSpace(1, datasetDims, maxdims);
        h5::DataSet *dataset = new h5::DataSet(group->createDataSet(setName, dtype, *dataSpace, prop));

        dataset->close();
        dataSpace->close();
        delete dataSpace;
        delete dataset;
    }
}

hsize_t utils::sizeOfFingerPrintDatasetInH5(h5::Group *group, const std::string setName) {
    if (!group->nameExists(setName)) {
        return 0;
    }
    
    hsize_t size[1];
    h5::DataSet *dataset = new h5::DataSet(group->openDataSet(setName));
    h5::DataSpace *fileSpace = new h5::DataSpace(dataset->getSpace());
    fileSpace->getSimpleExtentDims(size, NULL);

    dataset->close();
    fileSpace->close();
    delete fileSpace;
    delete dataset;
    return size[0];
}

bool utils::getFingerprintFromIndex(h5::Group *group, const std::string setName, hsize_t index, hsize_t CFPSize, uint64_t *cfp) {
    if (!group->nameExists(setName)) {
        return false;
    }

    // offset and size
    hsize_t offset[1] = {index};
    hsize_t datasetDims[1] = {CFPSize};

    //open dataset
    h5::DataSet *dataset = new h5::DataSet(group->openDataSet(setName));

    // open space in file
    h5::DataSpace *fileSpace = new h5::DataSpace(dataset->getSpace());

    // select hyperslab space in file
    fileSpace->selectHyperslab(H5S_SELECT_SET, datasetDims, offset);

    // read data
    h5::DataSpace *memSpace = new h5::DataSpace(1, datasetDims, NULL);
    dataset->read(cfp, dataset->getDataType(), *memSpace, *fileSpace);
    
    fileSpace->close();
    memSpace->close();
    dataset->close();
    delete fileSpace;
    delete memSpace;
    delete dataset;

    return true;
}

bool utils::updateFingerprintAtIndex(h5::Group *group, const std::string setName, hsize_t index, hsize_t CFPSize, uint64_t *cfp) {
    if (!group->nameExists(setName)) {
        return false;
    }

    // offset and size
    hsize_t offset[1] = {index};
    hsize_t datasetDims[1] = {CFPSize};

    //open dataset
    h5::DataSet *dataset = new h5::DataSet(group->openDataSet(setName));

    // open space in file and select hyperslab in file    
    h5::DataSpace *fileSpace = new h5::DataSpace(dataset->getSpace());
    fileSpace->selectHyperslab(H5S_SELECT_SET, datasetDims, offset);
    
    // write the modified data
    h5::DataSpace *memSpace = new h5::DataSpace(1, datasetDims, NULL);
    dataset->write(cfp, h5::PredType::NATIVE_UINT64, *memSpace, *fileSpace);

    fileSpace->close();
    memSpace->close();
    dataset->close();
    delete fileSpace;
    delete memSpace;
    delete dataset;
    return true;
}