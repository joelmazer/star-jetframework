#ifndef __STJETPICODEFINITIONS__H
#define __STJETPICODEFINITIONS__H

#include <TObject.h>
#include <iostream>

#define _NUMBER_OF_EVENT_FEATURES 7

class StJetPicoDefinitions : public TObject
{
 public:

  enum EVENT_FEATURES
  {
    _EVENT_WEIGHT          = 0, // this is for x-section weighting
    
    _EVENT_MEDIAN_PT       = 1, // from fj for example
    _EVENT_MEDIAN_PT_4VECT = 2,
    
    _EVENT_SIGMA           = 3,
    _EVENT_SIGMA_4VECT     = 4,
    
    _EVENT_MEAN_AREA       = 5,
    _EVENT_MEAN_AREA_4VECT = 6

  };
  
  static Int_t DebugLevel() {return gDebugLevel;}
  static void SetDebugLevel(Int_t level) 
    {
      gDebugLevel = level;
    }

 private:
  StJetPicoDefinitions() : TObject()
    {
      ;
    }

  virtual ~StJetPicoDefinitions()
    {
      ;
    }
  
  static Int_t gDebugLevel;
    
  ClassDef(StJetPicoDefinitions, 0)
};

#ifndef PICONODEBUG
#define __DEBUG(level, message) {if (level <= StJetPicoDefinitions::DebugLevel()) { std::cout << "[d] " << "[" << level << "] " << "[" << ClassName() << "::" << __FUNCTION__ << "] " << this->GetTitle() << ": " << message << std::endl;}}

#define __INFO(message) {std::cout << "[i] " << "[" << ClassName() << "::" << __FUNCTION__ << "] " << this->GetTitle() << ": " << message << std::endl;}

#define __ERROR(message) {std::cerr << "[e] " << "[" << ClassName() << "::" << __FUNCTION__ << "] " << this->GetTitle() << ": " << message << std::endl;}

#define __WARNING(message) {std::cout << "[w] " << "[" << ClassName() << "::" << __FUNCTION__ << "] " << this->GetTitle() << ": " << message << std::endl;}

#define __DEBUG_FUNCTION(level, message) {if (level <= StJetPicoDefinitions::DebugLevel()) { std::cout << "[d] " << "[" << level << "] " << "[" << __FUNCTION__ << "] " << message << std::endl;}}

#define __INFO_FUNCTION(message) {std::cout << "[i] [" << __FUNCTION__ << "] " << message << std::endl;}

#define __ERROR_FUNCTION(message) {std::cerr << "[e] [" << __FUNCTION__ << "] " << message << std::endl;}

#define __WARNING_FUNCTION(message) {std::cout << "[w] [" << __FUNCTION__ << "] " << message << std::endl;}
#endif

#ifdef PICONODEBUG
#define __DEBUG(level, message)
#define __INFO(message)
#define __ERROR(message)
#define __WARNING(message)
#define __DEBUG_FUNCTION(level, message)
#define __INFO_FUNCTION(message)
#define __ERROR_FUNCTION(message)
#define __WARNING_FUNCTION(message)
#endif

#endif
