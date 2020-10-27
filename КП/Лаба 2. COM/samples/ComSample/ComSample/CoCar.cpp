#include "CoCar.h"

// Конструктор и деструктор СоСаr

CoCar::CoCar()
	: m_refCount(0), m_currSpeed(0), m_maxSpeed(0)
{
	m_petName = SysAllocString(L"Default Pet Name");
}

CoCar::~CoCar()
{
	if (m_petName)
		SysFreeString(m_petName);
	MessageBox(NULL, "CoCar is dead", "Destructor", MB_OK |
		MB_SETFOREGROUND);
}

// IUnknown

STDMETHODIMP_(HRESULT __stdcall) CoCar::QueryInterface(REFIID riid, void** pIFace) {
	// Which aspect of me do they want?
	if (riid == IID_IUnknown)
	{
		*pIFace = (IUnknown*)(IEngine*)this;
		MessageBox(NULL, "Handed out IUnknown", "QI", MB_OK |
			MB_SETFOREGROUND);
	}

	else if (riid == IID_IEngine)
	{
		*pIFace = (IEngine*)this;
		MessageBox(NULL, "Handed out IEngine", "QI", MB_OK |
			MB_SETFOREGROUND);
	}

	else if (riid == IID_IStats)
	{
		*pIFace = (IStats*)this;
		MessageBox(NULL, "Handed out IStats", "QI", MB_OK |
			MB_SETFOREGROUND);
	}

	else if (riid == IID_ICreateCar)
	{
		*pIFace = (ICreateCar*)this;
		MessageBox(NULL, "Handed out ICreateCar", "QI", MB_OK |
			MB_SETFOREGROUND);
	}
	else
	{
		*pIFace = NULL;
		return E_NOINTERFACE;
	}

	((IUnknown*)(*pIFace))->AddRef();
	return S_OK;

}

STDMETHODIMP_(DWORD __stdcall) CoCar::AddRef() {
	++m_refCount;
	return m_refCount;
}

STDMETHODIMP_(DWORD __stdcall) CoCar::Release() {
	if (--m_refCount == 0)
	{
		delete this;
		return 0;
	}
	else
		return m_refCount;
}

// IEngine
// Реализация IEngine

STDMETHODIMP_(HRESULT __stdcall) CoCar::SpeedUp()
{
	m_currSpeed += 10;
	return S_OK;
}

STDMETHODIMP_(HRESULT __stdcall) CoCar::GetMaxSpeed(int* maxSpeed)
{
	*maxSpeed = m_maxSpeed;
	return S_OK;
}

STDMETHODIMP_(HRESULT __stdcall) CoCar::GetCurSpeed(int* curSpeed)
{
	*curSpeed = m_currSpeed;
	return S_OK;
}

// IStats
// Информация о СоСаr помещается в блоки сообщений

STDMETHODIMP_(HRESULT __stdcall) CoCar::DisplayStats()
{
	// Need to transfer a BSTR to a char array.
	char buff[MAX_NAME_LENGTH];
	WideCharToMultiByte(CP_ACP, NULL, m_petName, -1, buff,
		MAX_NAME_LENGTH, NULL, NULL);

	MessageBox(NULL, buff, "Pet Name", MB_OK | MB_SETFOREGROUND);
	memset(buff, 0, sizeof(buff));
	sprintf_s(buff, "%d", m_maxSpeed);
	MessageBox(NULL, buff, "Max Speed", MB_OK |
		MB_SETFOREGROUND);
	return S_OK;
}

// Возвращает клиенту копию внутреннего буфера BSTR 

STDMETHODIMP_(HRESULT __stdcall) CoCar::GetPetName(BSTR* petName)
{
	*petName = SysAllocString(m_petName);
	return S_OK;
}

// ICreateCar
// Реализация ICreateCar

STDMETHODIMP_(HRESULT __stdcall) CoCar::SetPetName(BSTR petName)
{
	SysReAllocString(&m_petName, petName);
	return S_OK;
}

STDMETHODIMP_(HRESULT __stdcall) CoCar::SetMaxSpeed(int maxSp)
{
	if (maxSp < MAX_SPEED)
		m_maxSpeed = maxSp;
	return S_OK;
}
