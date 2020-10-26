// Car.cpp : Implementation of CCar

#include "pch.h"
#include "Car.h"

// CCar
HRESULT CCar::FinalConstruct()
{
	_petName = SysAllocString(L"Default Pet Name");
	return S_OK;
}

void CCar::FinalRelease()
{
	if (_petName)
	{
		SysFreeString(_petName);
	}
}

STDMETHODIMP CCar::SpeedUp()
{
	_curSpeed += 10;
	return S_OK;
}

STDMETHODIMP CCar::GetMaxSpeed(int* maxSpeed)
{
	*maxSpeed = _maxSpeed;
	return S_OK;
}

STDMETHODIMP CCar::GetCurSpeed(int* curSpeed)
{
	*curSpeed = _curSpeed;
	return S_OK;
}

HRESULT CCar::DisplayStats()
{
	char buff[MAX_NAME_LENGTH];
	WideCharToMultiByte(CP_ACP, 0, _petName, -1, buff, MAX_NAME_LENGTH, nullptr, nullptr);
	MessageBoxA(nullptr, buff, "Pet name", MB_OK | MB_SETFOREGROUND);
	memset(buff, 0, sizeof buff);

	sprintf_s(buff, "%d", _maxSpeed);
	MessageBoxA(nullptr, buff, "Max speed", MB_OK | MB_SETFOREGROUND);

	return S_OK;
}

STDMETHODIMP CCar::GetPetName(BSTR* petName)
{
	*petName = SysAllocString(_petName);
	return S_OK;
}

STDMETHODIMP CCar::SetPetName(BSTR petName)
{
	SysReAllocString(&_petName, petName);
	return S_OK;
}

STDMETHODIMP CCar::SetMaxSpeed(int maxSp)
{
	if (maxSp < MAX_SPEED)
	{
		_maxSpeed = maxSp;
		return S_OK;
	}

	return E_INVALIDARG;
}
