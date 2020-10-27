#include "CoCarClassFactory.h"

ULONG g_lockCount = 0; // Количество блокировок сервера
ULONG g_objCount = 0; // Количество "живых" объектов в сервере


STDAPI DllGetClassObject(REFCLSID rclsid, REFIID riid, void** ppv)
{
	HRESULT hr;
	CoCarClassFactory* pCFact = NULL;

	if (rclsid != CLSID_CoCar)
		return CLASS_E_CLASSNOTAVAILABLE;

	pCFact = new CoCarClassFactory;

	MessageBox(NULL, "DllGetClassObject", "QI", MB_OK |
		MB_SETFOREGROUND);
	hr = pCFact->QueryInterface(riid, ppv);

	if (FAILED(hr))
		delete pCFact;

	return hr;
}

STDAPI DllCanUnloadNow()
{
	if (g_lockCount == 0 && g_objCount == 0)
	{
		return S_OK;
	}
	else
		return S_FALSE;
}
