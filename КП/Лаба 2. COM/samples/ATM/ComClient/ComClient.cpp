#include <iostream>
#include "../ComSample/interfaces.h"
#include "../ComSample/iid.h"

using namespace std;

int main()
{
	HRESULT hr;
	IClassFactory* pCF = NULL;
	ICreateAtm* pICreateAtm = NULL;
	IStats* pIStats = NULL;
	IBalance* pIBalance = NULL;

	hr = CoInitialize(NULL);
	if (!SUCCEEDED(hr))
		throw std::exception();

	hr = CoGetClassObject(CLSID_CoAtm, CLSCTX_INPROC_SERVER,
		NULL, IID_IClassFactory, (void**)&pCF);
	if (!SUCCEEDED(hr))
		throw std::exception();

	hr = pCF->CreateInstance(NULL, IID_ICreateAtm,
		(void**)&pICreateAtm);
	pCF->Release();
	if (!SUCCEEDED(hr))
		throw std::exception();

	if (SUCCEEDED(hr))
	{
		pICreateAtm->SetMoney(30000);
		BSTR addr = SysAllocString(L"Пермь, ул. Культуры, д. 10");
		pICreateAtm->SetLocationAddr(addr);
		SysFreeString(addr);

		// Now get IStats
		hr = pICreateAtm->QueryInterface(IID_IStats,
			(void**)&pIStats);
		pICreateAtm->Release();
	}

	if (SUCCEEDED(hr))
	{
		hr = pIStats->QueryInterface(IID_IBalance,
			(void**)&pIBalance);
	}

	if (SUCCEEDED(hr))
	{
		pIStats->DisplayStats();
		pIBalance->PayMoney(10252, 1000);
		pIBalance->PayMoney(10452, 3000);
		pIStats->DisplayStats();
		pIBalance->WithdrawMoney(10252, 500);
		pIStats->DisplayStats();

		if (pIBalance) pIBalance->Release();
		if (pIStats) pIStats->Release();
	}

	CoUninitialize();
	return 0;
}
