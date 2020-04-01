// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME myDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "Utils/TH1EFT.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *WCPoint_Dictionary();
   static void WCPoint_TClassManip(TClass*);
   static void *new_WCPoint(void *p = 0);
   static void *newArray_WCPoint(Long_t size, void *p);
   static void delete_WCPoint(void *p);
   static void deleteArray_WCPoint(void *p);
   static void destruct_WCPoint(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCPoint*)
   {
      ::WCPoint *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::WCPoint));
      static ::ROOT::TGenericClassInfo 
         instance("WCPoint", "Utils/WCPoint.h", 15,
                  typeid(::WCPoint), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &WCPoint_Dictionary, isa_proxy, 4,
                  sizeof(::WCPoint) );
      instance.SetNew(&new_WCPoint);
      instance.SetNewArray(&newArray_WCPoint);
      instance.SetDelete(&delete_WCPoint);
      instance.SetDeleteArray(&deleteArray_WCPoint);
      instance.SetDestructor(&destruct_WCPoint);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCPoint*)
   {
      return GenerateInitInstanceLocal((::WCPoint*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::WCPoint*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *WCPoint_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::WCPoint*)0x0)->GetClass();
      WCPoint_TClassManip(theClass);
   return theClass;
   }

   static void WCPoint_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *WCFit_Dictionary();
   static void WCFit_TClassManip(TClass*);
   static void *new_WCFit(void *p = 0);
   static void *newArray_WCFit(Long_t size, void *p);
   static void delete_WCFit(void *p);
   static void deleteArray_WCFit(void *p);
   static void destruct_WCFit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCFit*)
   {
      ::WCFit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::WCFit));
      static ::ROOT::TGenericClassInfo 
         instance("WCFit", "Utils/WCFit.h", 21,
                  typeid(::WCFit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &WCFit_Dictionary, isa_proxy, 4,
                  sizeof(::WCFit) );
      instance.SetNew(&new_WCFit);
      instance.SetNewArray(&newArray_WCFit);
      instance.SetDelete(&delete_WCFit);
      instance.SetDeleteArray(&deleteArray_WCFit);
      instance.SetDestructor(&destruct_WCFit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCFit*)
   {
      return GenerateInitInstanceLocal((::WCFit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::WCFit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *WCFit_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::WCFit*)0x0)->GetClass();
      WCFit_TClassManip(theClass);
   return theClass;
   }

   static void WCFit_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_TH1EFT(void *p = 0);
   static void *newArray_TH1EFT(Long_t size, void *p);
   static void delete_TH1EFT(void *p);
   static void deleteArray_TH1EFT(void *p);
   static void destruct_TH1EFT(void *p);
   static void directoryAutoAdd_TH1EFT(void *obj, TDirectory *dir);
   static Long64_t merge_TH1EFT(void *obj, TCollection *coll,TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TH1EFT*)
   {
      ::TH1EFT *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TH1EFT >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TH1EFT", ::TH1EFT::Class_Version(), "Utils/TH1EFT.h", 14,
                  typeid(::TH1EFT), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TH1EFT::Dictionary, isa_proxy, 4,
                  sizeof(::TH1EFT) );
      instance.SetNew(&new_TH1EFT);
      instance.SetNewArray(&newArray_TH1EFT);
      instance.SetDelete(&delete_TH1EFT);
      instance.SetDeleteArray(&deleteArray_TH1EFT);
      instance.SetDestructor(&destruct_TH1EFT);
      instance.SetDirectoryAutoAdd(&directoryAutoAdd_TH1EFT);
      instance.SetMerge(&merge_TH1EFT);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TH1EFT*)
   {
      return GenerateInitInstanceLocal((::TH1EFT*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TH1EFT*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TH1EFT::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TH1EFT::Class_Name()
{
   return "TH1EFT";
}

//______________________________________________________________________________
const char *TH1EFT::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TH1EFT*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TH1EFT::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TH1EFT*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TH1EFT::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TH1EFT*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TH1EFT::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TH1EFT*)0x0)->GetClass(); }
   return fgIsA;
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCPoint(void *p) {
      return  p ? new(p) ::WCPoint : new ::WCPoint;
   }
   static void *newArray_WCPoint(Long_t nElements, void *p) {
      return p ? new(p) ::WCPoint[nElements] : new ::WCPoint[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCPoint(void *p) {
      delete ((::WCPoint*)p);
   }
   static void deleteArray_WCPoint(void *p) {
      delete [] ((::WCPoint*)p);
   }
   static void destruct_WCPoint(void *p) {
      typedef ::WCPoint current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCPoint

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCFit(void *p) {
      return  p ? new(p) ::WCFit : new ::WCFit;
   }
   static void *newArray_WCFit(Long_t nElements, void *p) {
      return p ? new(p) ::WCFit[nElements] : new ::WCFit[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCFit(void *p) {
      delete ((::WCFit*)p);
   }
   static void deleteArray_WCFit(void *p) {
      delete [] ((::WCFit*)p);
   }
   static void destruct_WCFit(void *p) {
      typedef ::WCFit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCFit

//______________________________________________________________________________
void TH1EFT::Streamer(TBuffer &R__b)
{
   // Stream an object of class TH1EFT.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TH1EFT::Class(),this);
   } else {
      R__b.WriteClassBuffer(TH1EFT::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TH1EFT(void *p) {
      return  p ? new(p) ::TH1EFT : new ::TH1EFT;
   }
   static void *newArray_TH1EFT(Long_t nElements, void *p) {
      return p ? new(p) ::TH1EFT[nElements] : new ::TH1EFT[nElements];
   }
   // Wrapper around operator delete
   static void delete_TH1EFT(void *p) {
      delete ((::TH1EFT*)p);
   }
   static void deleteArray_TH1EFT(void *p) {
      delete [] ((::TH1EFT*)p);
   }
   static void destruct_TH1EFT(void *p) {
      typedef ::TH1EFT current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around the directory auto add.
   static void directoryAutoAdd_TH1EFT(void *p, TDirectory *dir) {
      ((::TH1EFT*)p)->DirectoryAutoAdd(dir);
   }
   // Wrapper around the merge function.
   static Long64_t  merge_TH1EFT(void *obj,TCollection *coll,TFileMergeInfo *) {
      return ((::TH1EFT*)obj)->Merge(coll);
   }
} // end of namespace ROOT for class ::TH1EFT

namespace ROOT {
   static TClass *vectorlEstringgR_Dictionary();
   static void vectorlEstringgR_TClassManip(TClass*);
   static void *new_vectorlEstringgR(void *p = 0);
   static void *newArray_vectorlEstringgR(Long_t size, void *p);
   static void delete_vectorlEstringgR(void *p);
   static void deleteArray_vectorlEstringgR(void *p);
   static void destruct_vectorlEstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<string>*)
   {
      vector<string> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<string>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<string>", -2, "vector", 216,
                  typeid(vector<string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEstringgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<string>) );
      instance.SetNew(&new_vectorlEstringgR);
      instance.SetNewArray(&newArray_vectorlEstringgR);
      instance.SetDelete(&delete_vectorlEstringgR);
      instance.SetDeleteArray(&deleteArray_vectorlEstringgR);
      instance.SetDestructor(&destruct_vectorlEstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<string> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<string>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<string>*)0x0)->GetClass();
      vectorlEstringgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEstringgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string> : new vector<string>;
   }
   static void *newArray_vectorlEstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string>[nElements] : new vector<string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEstringgR(void *p) {
      delete ((vector<string>*)p);
   }
   static void deleteArray_vectorlEstringgR(void *p) {
      delete [] ((vector<string>*)p);
   }
   static void destruct_vectorlEstringgR(void *p) {
      typedef vector<string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<string>

namespace ROOT {
   static TClass *vectorlEpairlEintcOintgRsPgR_Dictionary();
   static void vectorlEpairlEintcOintgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEpairlEintcOintgRsPgR(void *p = 0);
   static void *newArray_vectorlEpairlEintcOintgRsPgR(Long_t size, void *p);
   static void delete_vectorlEpairlEintcOintgRsPgR(void *p);
   static void deleteArray_vectorlEpairlEintcOintgRsPgR(void *p);
   static void destruct_vectorlEpairlEintcOintgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<pair<int,int> >*)
   {
      vector<pair<int,int> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<pair<int,int> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<pair<int,int> >", -2, "vector", 216,
                  typeid(vector<pair<int,int> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEpairlEintcOintgRsPgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<pair<int,int> >) );
      instance.SetNew(&new_vectorlEpairlEintcOintgRsPgR);
      instance.SetNewArray(&newArray_vectorlEpairlEintcOintgRsPgR);
      instance.SetDelete(&delete_vectorlEpairlEintcOintgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEpairlEintcOintgRsPgR);
      instance.SetDestructor(&destruct_vectorlEpairlEintcOintgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<pair<int,int> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<pair<int,int> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEpairlEintcOintgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<pair<int,int> >*)0x0)->GetClass();
      vectorlEpairlEintcOintgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEpairlEintcOintgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEpairlEintcOintgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<pair<int,int> > : new vector<pair<int,int> >;
   }
   static void *newArray_vectorlEpairlEintcOintgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<pair<int,int> >[nElements] : new vector<pair<int,int> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEpairlEintcOintgRsPgR(void *p) {
      delete ((vector<pair<int,int> >*)p);
   }
   static void deleteArray_vectorlEpairlEintcOintgRsPgR(void *p) {
      delete [] ((vector<pair<int,int> >*)p);
   }
   static void destruct_vectorlEpairlEintcOintgRsPgR(void *p) {
      typedef vector<pair<int,int> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<pair<int,int> >

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 216,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlEWCPointgR_Dictionary();
   static void vectorlEWCPointgR_TClassManip(TClass*);
   static void *new_vectorlEWCPointgR(void *p = 0);
   static void *newArray_vectorlEWCPointgR(Long_t size, void *p);
   static void delete_vectorlEWCPointgR(void *p);
   static void deleteArray_vectorlEWCPointgR(void *p);
   static void destruct_vectorlEWCPointgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<WCPoint>*)
   {
      vector<WCPoint> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<WCPoint>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<WCPoint>", -2, "vector", 216,
                  typeid(vector<WCPoint>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEWCPointgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<WCPoint>) );
      instance.SetNew(&new_vectorlEWCPointgR);
      instance.SetNewArray(&newArray_vectorlEWCPointgR);
      instance.SetDelete(&delete_vectorlEWCPointgR);
      instance.SetDeleteArray(&deleteArray_vectorlEWCPointgR);
      instance.SetDestructor(&destruct_vectorlEWCPointgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<WCPoint> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<WCPoint>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEWCPointgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<WCPoint>*)0x0)->GetClass();
      vectorlEWCPointgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEWCPointgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEWCPointgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<WCPoint> : new vector<WCPoint>;
   }
   static void *newArray_vectorlEWCPointgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<WCPoint>[nElements] : new vector<WCPoint>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEWCPointgR(void *p) {
      delete ((vector<WCPoint>*)p);
   }
   static void deleteArray_vectorlEWCPointgR(void *p) {
      delete [] ((vector<WCPoint>*)p);
   }
   static void destruct_vectorlEWCPointgR(void *p) {
      typedef vector<WCPoint> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<WCPoint>

namespace ROOT {
   static TClass *vectorlEWCFitgR_Dictionary();
   static void vectorlEWCFitgR_TClassManip(TClass*);
   static void *new_vectorlEWCFitgR(void *p = 0);
   static void *newArray_vectorlEWCFitgR(Long_t size, void *p);
   static void delete_vectorlEWCFitgR(void *p);
   static void deleteArray_vectorlEWCFitgR(void *p);
   static void destruct_vectorlEWCFitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<WCFit>*)
   {
      vector<WCFit> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<WCFit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<WCFit>", -2, "vector", 216,
                  typeid(vector<WCFit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEWCFitgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<WCFit>) );
      instance.SetNew(&new_vectorlEWCFitgR);
      instance.SetNewArray(&newArray_vectorlEWCFitgR);
      instance.SetDelete(&delete_vectorlEWCFitgR);
      instance.SetDeleteArray(&deleteArray_vectorlEWCFitgR);
      instance.SetDestructor(&destruct_vectorlEWCFitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<WCFit> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<WCFit>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEWCFitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<WCFit>*)0x0)->GetClass();
      vectorlEWCFitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEWCFitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEWCFitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<WCFit> : new vector<WCFit>;
   }
   static void *newArray_vectorlEWCFitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<WCFit>[nElements] : new vector<WCFit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEWCFitgR(void *p) {
      delete ((vector<WCFit>*)p);
   }
   static void deleteArray_vectorlEWCFitgR(void *p) {
      delete [] ((vector<WCFit>*)p);
   }
   static void destruct_vectorlEWCFitgR(void *p) {
      typedef vector<WCFit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<WCFit>

namespace ROOT {
   static TClass *unordered_maplEstringcOdoublegR_Dictionary();
   static void unordered_maplEstringcOdoublegR_TClassManip(TClass*);
   static void *new_unordered_maplEstringcOdoublegR(void *p = 0);
   static void *newArray_unordered_maplEstringcOdoublegR(Long_t size, void *p);
   static void delete_unordered_maplEstringcOdoublegR(void *p);
   static void deleteArray_unordered_maplEstringcOdoublegR(void *p);
   static void destruct_unordered_maplEstringcOdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const unordered_map<string,double>*)
   {
      unordered_map<string,double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(unordered_map<string,double>));
      static ::ROOT::TGenericClassInfo 
         instance("unordered_map<string,double>", -2, "unordered_map", 101,
                  typeid(unordered_map<string,double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &unordered_maplEstringcOdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(unordered_map<string,double>) );
      instance.SetNew(&new_unordered_maplEstringcOdoublegR);
      instance.SetNewArray(&newArray_unordered_maplEstringcOdoublegR);
      instance.SetDelete(&delete_unordered_maplEstringcOdoublegR);
      instance.SetDeleteArray(&deleteArray_unordered_maplEstringcOdoublegR);
      instance.SetDestructor(&destruct_unordered_maplEstringcOdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< unordered_map<string,double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const unordered_map<string,double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *unordered_maplEstringcOdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const unordered_map<string,double>*)0x0)->GetClass();
      unordered_maplEstringcOdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void unordered_maplEstringcOdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_unordered_maplEstringcOdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) unordered_map<string,double> : new unordered_map<string,double>;
   }
   static void *newArray_unordered_maplEstringcOdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) unordered_map<string,double>[nElements] : new unordered_map<string,double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_unordered_maplEstringcOdoublegR(void *p) {
      delete ((unordered_map<string,double>*)p);
   }
   static void deleteArray_unordered_maplEstringcOdoublegR(void *p) {
      delete [] ((unordered_map<string,double>*)p);
   }
   static void destruct_unordered_maplEstringcOdoublegR(void *p) {
      typedef unordered_map<string,double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class unordered_map<string,double>

namespace {
  void TriggerDictionaryInitialization_myDict_Impl() {
    static const char* headers[] = {
"Utils/TH1EFT.h",
0
    };
    static const char* includePaths[] = {
"/home/ntonon/Documents/Programmes/myBuild_root-6.18.04/include",
"/home/ntonon/Postdoc/TopEFT/Analysis/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "myDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$Utils/TH1EFT.h")))  WCPoint;
class __attribute__((annotate("$clingAutoload$Utils/TH1EFT.h")))  WCFit;
class __attribute__((annotate(R"ATTRDUMP(Needed to include custom class within ROOT)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$Utils/TH1EFT.h")))  TH1EFT;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "myDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "Utils/TH1EFT.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TH1EFT", payloadCode, "@",
"WCFit", payloadCode, "@",
"WCPoint", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("myDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_myDict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_myDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_myDict() {
  TriggerDictionaryInitialization_myDict_Impl();
}
