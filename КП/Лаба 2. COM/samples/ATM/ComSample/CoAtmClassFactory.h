#include "CoAtm.h"

class CoAtmClassFactory : public IClassFactory
{
public:
	CoAtmClassFactory();

	virtual ~CoAtmClassFactory();

	// IUnknown
	STDMETHODIMP_(ULONG) AddRef();

	STDMETHODIMP_(ULONG) Release();

	STDMETHODIMP QueryInterface(REFIID riid, void** ppv);


	// ICF
	STDMETHODIMP LockServer(BOOL fLock);

	STDMETHODIMP CreateInstance(LPUNKNOWN pUnkOuter, REFIID riid, void** ppv);


private:
	DWORD m_refCount;
};

