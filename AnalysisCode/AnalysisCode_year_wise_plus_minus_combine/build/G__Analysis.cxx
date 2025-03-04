// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__Analysis
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
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

// Header files passed as explicit arguments
#include "/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/AnalysisCode/AnalysisCode_year_wise_plus_minus_combine/CumulativePlots.h"
#include "/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/AnalysisCode/AnalysisCode_year_wise_plus_minus_combine/individual_channels.h"
#include "/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/AnalysisCode/AnalysisCode_year_wise_plus_minus_combine/MeanPlots.h"
#include "/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/AnalysisCode/AnalysisCode_year_wise_plus_minus_combine/SystematicRemoval.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *individual_channels_Dictionary();
   static void individual_channels_TClassManip(TClass*);
   static void *new_individual_channels(void *p = nullptr);
   static void *newArray_individual_channels(Long_t size, void *p);
   static void delete_individual_channels(void *p);
   static void deleteArray_individual_channels(void *p);
   static void destruct_individual_channels(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::individual_channels*)
   {
      ::individual_channels *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::individual_channels));
      static ::ROOT::TGenericClassInfo 
         instance("individual_channels", "individual_channels.h", 12,
                  typeid(::individual_channels), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &individual_channels_Dictionary, isa_proxy, 0,
                  sizeof(::individual_channels) );
      instance.SetNew(&new_individual_channels);
      instance.SetNewArray(&newArray_individual_channels);
      instance.SetDelete(&delete_individual_channels);
      instance.SetDeleteArray(&deleteArray_individual_channels);
      instance.SetDestructor(&destruct_individual_channels);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::individual_channels*)
   {
      return GenerateInitInstanceLocal(static_cast<::individual_channels*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::individual_channels*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *individual_channels_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::individual_channels*>(nullptr))->GetClass();
      individual_channels_TClassManip(theClass);
   return theClass;
   }

   static void individual_channels_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *SystematicRemoval_Dictionary();
   static void SystematicRemoval_TClassManip(TClass*);
   static void *new_SystematicRemoval(void *p = nullptr);
   static void *newArray_SystematicRemoval(Long_t size, void *p);
   static void delete_SystematicRemoval(void *p);
   static void deleteArray_SystematicRemoval(void *p);
   static void destruct_SystematicRemoval(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SystematicRemoval*)
   {
      ::SystematicRemoval *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::SystematicRemoval));
      static ::ROOT::TGenericClassInfo 
         instance("SystematicRemoval", "SystematicRemoval.h", 12,
                  typeid(::SystematicRemoval), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &SystematicRemoval_Dictionary, isa_proxy, 0,
                  sizeof(::SystematicRemoval) );
      instance.SetNew(&new_SystematicRemoval);
      instance.SetNewArray(&newArray_SystematicRemoval);
      instance.SetDelete(&delete_SystematicRemoval);
      instance.SetDeleteArray(&deleteArray_SystematicRemoval);
      instance.SetDestructor(&destruct_SystematicRemoval);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SystematicRemoval*)
   {
      return GenerateInitInstanceLocal(static_cast<::SystematicRemoval*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::SystematicRemoval*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *SystematicRemoval_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::SystematicRemoval*>(nullptr))->GetClass();
      SystematicRemoval_TClassManip(theClass);
   return theClass;
   }

   static void SystematicRemoval_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_individual_channels(void *p) {
      return  p ? new(p) ::individual_channels : new ::individual_channels;
   }
   static void *newArray_individual_channels(Long_t nElements, void *p) {
      return p ? new(p) ::individual_channels[nElements] : new ::individual_channels[nElements];
   }
   // Wrapper around operator delete
   static void delete_individual_channels(void *p) {
      delete (static_cast<::individual_channels*>(p));
   }
   static void deleteArray_individual_channels(void *p) {
      delete [] (static_cast<::individual_channels*>(p));
   }
   static void destruct_individual_channels(void *p) {
      typedef ::individual_channels current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::individual_channels

namespace ROOT {
   // Wrappers around operator new
   static void *new_SystematicRemoval(void *p) {
      return  p ? new(p) ::SystematicRemoval : new ::SystematicRemoval;
   }
   static void *newArray_SystematicRemoval(Long_t nElements, void *p) {
      return p ? new(p) ::SystematicRemoval[nElements] : new ::SystematicRemoval[nElements];
   }
   // Wrapper around operator delete
   static void delete_SystematicRemoval(void *p) {
      delete (static_cast<::SystematicRemoval*>(p));
   }
   static void deleteArray_SystematicRemoval(void *p) {
      delete [] (static_cast<::SystematicRemoval*>(p));
   }
   static void destruct_SystematicRemoval(void *p) {
      typedef ::SystematicRemoval current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::SystematicRemoval

namespace {
  void TriggerDictionaryInitialization_libAnalysis_Impl() {
    static const char* headers[] = {
"/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/AnalysisCode/AnalysisCode_year_wise_plus_minus_combine/CumulativePlots.h",
"/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/AnalysisCode/AnalysisCode_year_wise_plus_minus_combine/individual_channels.h",
"/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/AnalysisCode/AnalysisCode_year_wise_plus_minus_combine/MeanPlots.h",
"/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/AnalysisCode/AnalysisCode_year_wise_plus_minus_combine/SystematicRemoval.h",
nullptr
    };
    static const char* includePaths[] = {
"/usr/include/root",
"/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/AnalysisCode/AnalysisCode_year_wise_plus_minus_combine",
"/usr/include/root",
"/afs/cern.ch/work/n/nrawal/CSCAgeing_code_study/AnalysisCode/AnalysisCode_year_wise_plus_minus_combine/build/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libAnalysis dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/AnalysisCode/AnalysisCode_year_wise_plus_minus_combine/individual_channels.h")))  individual_channels;
class __attribute__((annotate("$clingAutoload$/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/AnalysisCode/AnalysisCode_year_wise_plus_minus_combine/SystematicRemoval.h")))  SystematicRemoval;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libAnalysis dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/AnalysisCode/AnalysisCode_year_wise_plus_minus_combine/CumulativePlots.h"
#include "/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/AnalysisCode/AnalysisCode_year_wise_plus_minus_combine/individual_channels.h"
#include "/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/AnalysisCode/AnalysisCode_year_wise_plus_minus_combine/MeanPlots.h"
#include "/afs/cern.ch/user/n/nrawal/work/CSCAgeing_code_study/AnalysisCode/AnalysisCode_year_wise_plus_minus_combine/SystematicRemoval.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"SystematicRemoval", payloadCode, "@",
"individual_channels", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libAnalysis",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libAnalysis_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libAnalysis_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libAnalysis() {
  TriggerDictionaryInitialization_libAnalysis_Impl();
}
