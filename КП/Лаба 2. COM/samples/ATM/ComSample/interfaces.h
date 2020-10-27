#pragma once

#ifndef _INTERFACES
#define _INTERFACES

#include <windows.h>

//IStats используется для получения адреса банкомата и статистики по взятию/вводу средств
DECLARE_INTERFACE_(IStats, IUnknown)
{
	STDMETHOD(DisplayStats)() PURE;
	STDMETHOD(GetLocationAddr)(BSTR * addr) PURE;
};

//IBalance позволяет управлять балансом владельцев счетов
DECLARE_INTERFACE_(IBalance, IUnknown)
{
	STDMETHOD(PayMoney)(int accountId, int money) PURE;
	STDMETHOD(WithdrawMoney)(int accountId, int money) PURE;
	STDMETHOD(GetBalance)(int accountId, int* balance) PURE;
};

//ICreateAtm используется для присвоения адреса банкомату и установке начального кол-ва средств
DECLARE_INTERFACE_(ICreateAtm, IUnknown)
{
	STDMETHOD(SetLocationAddr)(BSTR addr) PURE;
	STDMETHOD(SetMoney)(int money) PURE;
};

#endif // _INTERFACES
