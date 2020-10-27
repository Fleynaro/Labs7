#include "CoAtmClassFactory.h"
#include <string>

ULONG g_lockCount = 0; // Количество блокировок сервера
ULONG g_objCount = 0; // Количество "живых" объектов в сервере

STDAPI DllGetClassObject(REFCLSID rclsid, REFIID riid, void** ppv)
{
	HRESULT hr;
	CoAtmClassFactory* pCFact = NULL;

	if (rclsid != CLSID_CoAtm)
		return CLASS_E_CLASSNOTAVAILABLE;

	pCFact = new CoAtmClassFactory;

	MessageBox(NULL, (std::string("DllGetClassObject, g_objCount = ") + std::to_string(g_objCount)).c_str(), "QI", MB_OK |
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
