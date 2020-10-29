#include "CoAtm.h"

extern DWORD g_objCount;

#include <string>

CoAtm::CoAtm()
	: m_refCount(0), m_money(0)
{
	m_locAddr = SysAllocString(L"Default location address");
	g_objCount++;
}

CoAtm::~CoAtm()
{
	--g_objCount;
	if (m_locAddr)
		SysFreeString(m_locAddr);
	/*MessageBox(NULL, "CoAtm is dead", "Destructor", MB_OK |
		MB_SETFOREGROUND);*/
}

// IUnknown

STDMETHODIMP_(HRESULT __stdcall) CoAtm::QueryInterface(REFIID riid, void** pIFace) {
	/*MessageBox(NULL, (std::string("QueryInterface = ") + std::to_string(riid.Data1)).c_str(), "QI", MB_OK |
	MB_SETFOREGROUND);*/

	// Which aspect of me do they want?
	if (riid == IID_IUnknown)
	{
		*pIFace = (IUnknown*)(IBalance*)this;
		/*MessageBox(NULL, "Handed out IUnknown", "QI", MB_OK |
			MB_SETFOREGROUND);*/
	}

	else if (riid == IID_IBalance)
	{
		*pIFace = (IBalance*)this;
		/*MessageBox(NULL, "Handed out IBalance", "QI", MB_OK |
			MB_SETFOREGROUND);*/
	}

	else if (riid == IID_IStats)
	{
		*pIFace = (IStats*)this;
		/*MessageBox(NULL, "Handed out IStats", "QI", MB_OK |
			MB_SETFOREGROUND);*/
	}

	else if (riid == IID_ICreateAtm)
	{
		*pIFace = (ICreateAtm*)this;
		/*MessageBox(NULL, "Handed out ICreateAtm", "QI", MB_OK |
			MB_SETFOREGROUND);*/
	}
	else
	{
		*pIFace = NULL;
		return E_NOINTERFACE;
	}

	((IUnknown*)(*pIFace))->AddRef();
	return S_OK;

}

STDMETHODIMP_(DWORD __stdcall) CoAtm::AddRef() {
	++m_refCount;
	return m_refCount;
}

STDMETHODIMP_(DWORD __stdcall) CoAtm::Release() {
	if (--m_refCount == 0)
	{
		delete this;
		return 0;
	}
	else
		return m_refCount;
}

// IStats
// Информация о СоСаr помещается в блоки сообщений

STDMETHODIMP_(HRESULT __stdcall) CoAtm::DisplayStats()
{
	// Need to transfer a BSTR to a char array.
	char buff[MAX_NAME_LENGTH];
	WideCharToMultiByte(CP_ACP, NULL, m_locAddr, -1, buff,
		MAX_NAME_LENGTH, NULL, NULL);

	MessageBox(NULL, buff, "Loc addr", MB_OK | MB_SETFOREGROUND);
	memset(buff, 0, sizeof(buff));
	sprintf_s(buff, "%d", m_money);
	MessageBox(NULL, buff, "Money", MB_OK |
		MB_SETFOREGROUND);
	return S_OK;
}

// Возвращает клиенту копию внутреннего буфера BSTR 

STDMETHODIMP_(HRESULT __stdcall) CoAtm::GetLocationAddr(BSTR* addr)
{
	*addr = SysAllocString(m_locAddr);
	return S_OK;
}

// ICreateAtm
// Реализация ICreateAtm

STDMETHODIMP_(HRESULT __stdcall) CoAtm::SetLocationAddr(BSTR addr)
{
	SysReAllocString(&m_locAddr, addr);
	return S_OK;
}

// Inherited via IBalance

STDMETHODIMP_(HRESULT __stdcall) CoAtm::PayMoney(int accountId, int money)
{
	m_money += money;
	return S_OK;
}

STDMETHODIMP_(HRESULT __stdcall) CoAtm::WithdrawMoney(int accountId, int money)
{
	if (money > m_money)
		return E_FAIL;
	m_money -= money;
	return S_OK;
}

STDMETHODIMP_(HRESULT __stdcall) CoAtm::GetBalance(int accountId, int* balance)
{
	*balance = 0;
	return S_OK;
}

// Inherited via ICreateAtm

STDMETHODIMP_(HRESULT __stdcall) CoAtm::SetMoney(int money)
{
	m_money = money;
	return S_OK;
}

STDMETHODIMP CoAtm::GetMoney(int* money)
{
	*money = m_money;
	return S_OK;
}
