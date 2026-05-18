#ifndef ABSTRACT_EVENT_SAVER
#define ABSTRACT_EVENT_SAVER

#include "Event.h"

class AbstractEventSaver {
  public:
    AbstractEventSaver();
    virtual ~AbstractEventSaver();
    virtual void Save(Event *, Model *, int) = 0;
    void FindPreviousEventFiles();
    void SetEventsTemp();

  protected:
    void ReadParameters();
    void SaveParameters(Model *tModel);
    int	        mFileCounter;
    const int   kEventsPerFile;
    int	        mNumberOfEvents;

    TTree*	        mParameterTree;

    Bool_t          mRandomize;
    Int_t           mMaxIntegrationSamples;
    Double_t        mIntegTolerance;
    Int_t           mIntegToleranceInterval;
    Int_t           mIntegToleranceNSuccessive;

};

#endif
