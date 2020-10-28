#pragma once

#ifndef _INTERFACES
#define _INTERFACES

#include <windows.h>

//IStats ������������ ��� ��������� ������ ��������� � ���������� �� ������/����� �������
DECLARE_INTERFACE_(IStats, IUnknown)
{
	STDMETHOD(DisplayStats)() PURE;
	STDMETHOD(GetLocationAddr)(BSTR * addr) PURE;
};

//IBalance ��������� ��������� �������� ���������� ������
DECLARE_INTERFACE_(IBalance, IUnknown)
{
	STDMETHOD(PayMoney)(int accountId, int money) PURE;
	STDMETHOD(WithdrawMoney)(int accountId, int money) PURE;
	STDMETHOD(GetBalance)(int accountId, int* balance) PURE;
};

//ICreateAtm ������������ ��� ���������� ������ ��������� � ��������� ���������� ���-�� �������
DECLARE_INTERFACE_(ICreateAtm, IUnknown)
{
	STDMETHOD(SetLocationAddr)(BSTR addr) PURE;
	STDMETHOD(SetMoney)(int money) PURE;
};

#endif // _INTERFACES
