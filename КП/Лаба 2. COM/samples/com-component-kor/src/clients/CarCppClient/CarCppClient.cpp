#include <iostream>
#include <Windows.h>
#include "../../components/CarComponent/CarComponent_i.h"
#include "../../components/CarComponent/CarComponent_i.c"

int main()
{
	CoInitialize(nullptr);

	IClassFactory* pCF = nullptr;
	CoGetClassObject(CLSID_Car, CLSCTX_INPROC_SERVER, nullptr, IID_IClassFactory,
		reinterpret_cast<void**>(&pCF));

	ICreateCar* pICreateCar = nullptr;
	auto hr = pCF->CreateInstance(nullptr, IID_ICreateCar, reinterpret_cast<void**>(&pICreateCar));

	pCF->Release();

	IStats* pStats = nullptr;
	if (SUCCEEDED(hr))
	{
		pICreateCar->SetMaxSpeed(30);
		auto* const petName = SysAllocString(L"Hedgehog John");
		pICreateCar->SetPetName(petName);
		SysFreeString(petName);

		hr = pICreateCar->QueryInterface(IID_IStats, reinterpret_cast<void**>(&pStats));
		pICreateCar->Release();
	}

	IEngine* pEngine = nullptr;
	if (SUCCEEDED(hr))
	{
		pStats->DisplayStats();
		hr = pStats->QueryInterface(IID_IEngine, reinterpret_cast<void**>(&pEngine));
	}

	if (SUCCEEDED(hr))
	{
		auto maxSp = 0;
		pEngine->GetMaxSpeed(&maxSp);

		auto curSp = 0;
		do
		{
			pEngine->SpeedUp();

			pEngine->GetCurSpeed(&curSp);
			std::cout << "Speed is " << curSp << " kph\n";
		}
		while (curSp <= maxSp);

		char buff[128];
		BSTR bstr;
		pStats->GetPetName(&bstr);
		WideCharToMultiByte(CP_ACP, 0, bstr, -1, buff, 128, nullptr, nullptr);

		std::cout << buff << " has blown up!\n\n";

		SysFreeString(bstr);

		if (pEngine)
		{
			pEngine->Release();
		}

		if (pStats)
		{
			pStats->Release();
		}
	}

	CoUninitialize();
}
