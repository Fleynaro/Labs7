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

	//class DispatchImpl : public IDispatch
	//{
	//public:
	//	// IDispatch methods
	//	STDMETHODIMP GetTypeInfoCount(
	//		PUINT			pctinfo)
	//	{
	//		*pctinfo = 1;
	//		return S_OK;
	//	}

	//	STDMETHODIMP GetTypeInfo(
	//		UINT			iTInfo,
	//		LCID			lcid,
	//		ITypeInfo** ppTInfo)
	//	{
	//		HRESULT					hr;
	//		ITypeLib* pTypeLib;

	//		*ppTInfo = NULL;
	//		if (0 != iTInfo) return DISP_E_BADINDEX;
	//		if (NULL == this->m_pTypeInfo)
	//		{
	//			hr = GetServer()->GetTypeLib(&pTypeLib);
	//			if (SUCCEEDED(hr))
	//			{
	//				hr = pTypeLib->GetTypeInfoOfGuid(MY_IID, &(this->m_pTypeInfo));
	//				pTypeLib->Release();
	//			}
	//		}
	//		else hr = S_OK;
	//		if (SUCCEEDED(hr))
	//		{
	//			*ppTInfo = this->m_pTypeInfo;
	//			this->m_pTypeInfo->AddRef();
	//		}
	//		return hr;
	//	}

	//	STDMETHODIMP GetIDsOfNames(
	//		REFIID			riid,
	//		OLECHAR** rgszNames,
	//		UINT			cNames,
	//		LCID			lcid,
	//		DISPID* rgDispId)
	//	{
	//		HRESULT					hr;
	//		ITypeInfo* pTypeInfo;
	//		hr = this->GetTypeInfo(0, LOCALE_SYSTEM_DEFAULT, &pTypeInfo);
	//		if (SUCCEEDED(hr))
	//		{
	//			hr = DispGetIDsOfNames(pTypeInfo, rgszNames, cNames, rgDispId);
	//			pTypeInfo->Release();
	//		}
	//		return hr;
	//	}

	//	STDMETHODIMP Invoke(
	//		DISPID			dispIdMember,
	//		REFIID			riid,
	//		LCID			lcid,
	//		WORD			wFlags,
	//		DISPPARAMS* pDispParams,
	//		VARIANT* pVarResult,
	//		EXCEPINFO* pExcepInfo,
	//		PUINT			puArgErr)
	//	{

	//	}

	//	HRESULT GetTestFunc(
	//		DISPPARAMS* pDispParams,
	//		VARIANT* pVarResult)
	//	{
	//		if (NULL == pVarResult) return E_INVALIDARG;
	//		HRESULT				hr;
	//		VARIANTARG			varg;
	//		UINT				uArgErr;
	//		VariantInit(&varg);
	//		hr = DispGetParam(pDispParams, 0, VT_I4, &varg, &uArgErr);
	//		if (FAILED(hr)) return hr;
	//		VariantInit(pVarResult);
	//		pVarResult->vt = VT_DISPATCH;
	//		
	//		return S_OK;
	//	}
	//};

};
