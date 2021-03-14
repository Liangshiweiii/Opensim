#ifndef OPOBJECT_H
#define OPOBJECT_H

#include "Tools/Includes.h"

namespace opensim
{

class Settings;

class OPObject
{
 public:
    virtual ~OPObject(void);
    std::string thisclassname;                                                  ///< Object's class name对象的类名
    std::string thisobjectname;                                                 ///< Object's name对象名
    //std::string DefaultInputFileName;											///< Default input file name默认输入文件名

    virtual void Initialize(const Settings& locSettings){};
    virtual void ReadInput(void);                                               ///< Reads input data from the default input file从默认输入文件中读入输入数据

    virtual void ReadInput(std::string InputFileName){};                        ///< Reads input data from the user specified input file. Should be implemented in each derived class从用户指定输入文件中读入输入数据，需要在每个衍生类中进行补充。
    bool IsInitialized(void) const;                                             ///< Returns True if Initialize() method has been called, False otherwise如果Initial（）方法被使用时，返回true，否则返回false
    bool IsNotInitialized(void) const;                                          ///< Returns False if Initialize() method has been called, True otherwise如果Initial（）方法被使用时返回false，否则返回true
    virtual void Remesh(size_t newNx, size_t newNy, size_t newNz){};            ///< Remeshes the object's storages重新划分对象的储存空间
    static OPObject* findOPObject(std::vector<OPObject*> ObjectList,
                           std::string ObjectName, std::string thisclassname,
                           std::string thisfunctionname, bool necessary, bool verbose = true)
    /* Returns a reference to the first OPObject in the list whose class name
     * starts with ObjectsName. I.e.: if ObjectsName is "Elasticity", the first
     * object in the list that has either the name ElasticitySteinbach or
     * ElasticityKhachaturyan will be returned.
     */
    {
        for(unsigned int i = 0; i < ObjectList.size(); i++)
        {
            if(ObjectList[i]->thisclassname.size() >= ObjectName.size())
            {
                std::string locObjectName = ObjectList[i]->thisclassname;
                locObjectName.resize(ObjectName.size());

                if(locObjectName == ObjectName)
                    return ObjectList[i];
            }
        }
        if (verbose)
        {
            std::cout << "No " << ObjectName << " object found! " << thisclassname
                      << "."<<thisfunctionname <<"()"<< std::endl;
        }

        if(necessary)
        {
            exit(13);
        }

        return nullptr;
    }

    static std::vector<OPObject*> findOPObjects(std::vector<OPObject*> ObjectList,
                           std::string ObjectName, std::string thisclassname,
                           std::string thisfunctionname, bool necessary, bool verbose = true)
    /* Returns a list with OPObject references found in the passed ObjectList
     * whose class name start with ObjectName. I.e.: if ObjectsName is
     * "Elasticity", all objects in the ObjectList that have a name starting
     * with "Elasticity", e.g. ElasticitySteinbach or ElasticityKhachaturyan
     * will be referred in the result list.
     */
    {
        std::vector<OPObject*> result;

        for(unsigned int i = 0; i < ObjectList.size(); i++)
        {
            if(ObjectList[i]->thisclassname.size() >= ObjectName.size())
            {
                std::string locObjectName = ObjectList[i]->thisclassname;
                locObjectName.resize(ObjectName.size());

                if(locObjectName == ObjectName)
                    result.push_back(ObjectList[i]);
            }
        }

        if(result.size()==0)
        {
            if (verbose)
            {
                std::cout << "No " << ObjectName << " object found! " << thisclassname
                          << thisfunctionname << std::endl;
            }
            if(necessary)
            {
                exit(13);
            }
        }

        return result;
    }

    bool initialized = false;                                                   ///< Object initialization status flag对象初始化状态标签
 protected:
 private:
};

}// namespace opensim
#endif