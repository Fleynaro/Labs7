#include "interfaces.h"
#include "iid.h"
#include <stdio.h>
#include <map>

#define MAX_NAME_LENGTH 64
#define MAX_SPEED 100

class CoAtm : public IBalance, public ICreateAtm, public IStats
{
	int m_refCount;
	BSTR	m_locAddr; // Инициализация через SysAllocString(), 
	// удаление — через вызов SysFreeString()
	int		m_money; // Деньги
public:
	// Конструктор и деструктор СоСаr
	CoAtm();

	virtual ~CoAtm();
	
	// IUnknown
	STDMETHODIMP QueryInterface(REFIID riid, void** pIFace) override;

	STDMETHODIMP_(DWORD)AddRef() override;

	STDMETHODIMP_(DWORD)Release() override;
	
	// IStats
	// Информация о СоСаr помещается в блоки сообщений
	STDMETHODIMP DisplayStats() override;

	// Возвращает клиенту копию внутреннего буфера BSTR 
	STDMETHODIMP GetLocationAddr(BSTR* addr) override;

	// ICreateAtm
	// Реализация ICreateAtm
	STDMETHODIMP SetLocationAddr(BSTR addr) override;


	// Inherited via IBalance
	STDMETHODIMP PayMoney(int accountId, int money) override;

	STDMETHODIMP WithdrawMoney(int accountId, int money) override;

	STDMETHODIMP GetBalance(int accountId, int* balance) override;

	// Inherited via ICreateAtm
	STDMETHODIMP SetMoney(int money) override;

};
