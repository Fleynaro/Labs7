#include "CoAtmClassFactory.h"

extern DWORD g_lockCount;
extern DWORD g_objCount;

CoAtmClassFactory::CoAtmClassFactory()
{
	m_refCount = 0;
	g_objCount++;
}

CoAtmClassFactory::~CoAtmClassFactory()
{
	g_objCount--;
}

// IUnknown

STDMETHODIMP_(ULONG __stdcall) CoAtmClassFactory::AddRef()
{
	return ++m_refCount;

}

STDMETHODIMP_(ULONG __stdcall) CoAtmClassFactory::Release()
{
	if (--m_refCount == 0)
	{
		delete this;
		return 0;
	}
	return m_refCount;
}

STDMETHODIMP_(HRESULT __stdcall) CoAtmClassFactory::QueryInterface(REFIID riid, void** ppv)
{
	// Which aspect of me do they want?
	if (riid == IID_IUnknown)
	{
		*ppv = (IUnknown*)this;
	}
	else if (riid == IID_IClassFactory)
	{
		*ppv = (IClassFactory*)this;
	}
	else
	{
		*ppv = NULL;
		return E_NOINTERFACE;
	}

	((IUnknown*)(*ppv))->AddRef();
	return S_OK;
}

// ICF

STDMETHODIMP_(HRESULT __stdcall) CoAtmClassFactory::LockServer(BOOL fLock)
{
	if (fLock)
		++g_lockCount;
	else
		--g_lockCount;

	return S_OK;
}

STDMETHODIMP_(HRESULT __stdcall) CoAtmClassFactory::CreateInstance(LPUNKNOWN pUnkOuter, REFIID riid, void** ppv)
{
	if (pUnkOuter != NULL)
		return CLASS_E_NOAGGREGATION;

	CoAtm* pCarObj = NULL;
	HRESULT hr;

	pCarObj = new CoAtm;

	hr = pCarObj->QueryInterface(riid, ppv);

	if (FAILED(hr))
		delete pCarObj;

	return hr;	// S_OK or E_NOINTERFACE.

}
