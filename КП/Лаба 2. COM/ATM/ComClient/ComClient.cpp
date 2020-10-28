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
		system("chcp 1251");
		printf("Введите номер карты:\n");
		int accountId;
		cin >> accountId;

		while (true) {
			printf("Введите:\n1, чтобы положить деньги\n2, чтобы снять деньги\n3, чтобы узнать, сколько наличных в банкомате");
			int actionId;
			cin >> actionId;

			if (actionId == 1 || actionId == 2) {
				printf("Введите сумму:\n");
				int money;
				cin >> money;
				if (actionId == 1) {
					pIBalance->PayMoney(accountId, money);
				}
				else {
					if (pIBalance->WithdrawMoney(accountId, money) == E_FAIL) {
						printf("Error: Не хватает денег в банкомате\n");
					}
				}
			}
			else if (actionId == 3) {
				pIStats->DisplayStats();
			}
		}

		
		if (pIBalance) pIBalance->Release();
		if (pIStats) pIStats->Release();
	}

	CoUninitialize();
	return 0;
}
