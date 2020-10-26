

/* this ALWAYS GENERATED file contains the definitions for the interfaces */


 /* File created by MIDL compiler version 8.01.0622 */
/* at Tue Jan 19 08:14:07 2038
 */
/* Compiler settings for CarComponent.idl:
    Oicf, W1, Zp8, env=Win32 (32b run), target_arch=X86 8.01.0622 
    protocol : dce , ms_ext, c_ext, robust
    error checks: allocation ref bounds_check enum stub_data 
    VC __declspec() decoration level: 
         __declspec(uuid()), __declspec(selectany), __declspec(novtable)
         DECLSPEC_UUID(), MIDL_INTERFACE()
*/
/* @@MIDL_FILE_HEADING(  ) */



/* verify that the <rpcndr.h> version is high enough to compile this file*/
#ifndef __REQUIRED_RPCNDR_H_VERSION__
#define __REQUIRED_RPCNDR_H_VERSION__ 500
#endif

#include "rpc.h"
#include "rpcndr.h"

#ifndef __RPCNDR_H_VERSION__
#error this stub requires an updated version of <rpcndr.h>
#endif /* __RPCNDR_H_VERSION__ */

#ifndef COM_NO_WINDOWS_H
#include "windows.h"
#include "ole2.h"
#endif /*COM_NO_WINDOWS_H*/

#ifndef __CarComponent_i_h__
#define __CarComponent_i_h__

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

/* Forward Declarations */ 

#ifndef __IEngine_FWD_DEFINED__
#define __IEngine_FWD_DEFINED__
typedef interface IEngine IEngine;

#endif 	/* __IEngine_FWD_DEFINED__ */


#ifndef __IStats_FWD_DEFINED__
#define __IStats_FWD_DEFINED__
typedef interface IStats IStats;

#endif 	/* __IStats_FWD_DEFINED__ */


#ifndef __ICreateCar_FWD_DEFINED__
#define __ICreateCar_FWD_DEFINED__
typedef interface ICreateCar ICreateCar;

#endif 	/* __ICreateCar_FWD_DEFINED__ */


#ifndef __ICar_FWD_DEFINED__
#define __ICar_FWD_DEFINED__
typedef interface ICar ICar;

#endif 	/* __ICar_FWD_DEFINED__ */


#ifndef __Car_FWD_DEFINED__
#define __Car_FWD_DEFINED__

#ifdef __cplusplus
typedef class Car Car;
#else
typedef struct Car Car;
#endif /* __cplusplus */

#endif 	/* __Car_FWD_DEFINED__ */


/* header files for imported files */
#include "oaidl.h"
#include "ocidl.h"
#include "shobjidl.h"

#ifdef __cplusplus
extern "C"{
#endif 


#ifndef __IEngine_INTERFACE_DEFINED__
#define __IEngine_INTERFACE_DEFINED__

/* interface IEngine */
/* [unique][helpstring][dual][uuid][object] */ 


EXTERN_C const IID IID_IEngine;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("9E37E773-86BD-4DD5-8F74-30536AAFC3B3")
    IEngine : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE SpeedUp( void) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetMaxSpeed( 
            /* [out] */ int *maxSpeed) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetCurSpeed( 
            /* [out] */ int *curSpeed) = 0;
        
    };
    
    
#else 	/* C style interface */

    typedef struct IEngineVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE *QueryInterface )( 
            IEngine * This,
            /* [in] */ REFIID riid,
            /* [annotation][iid_is][out] */ 
            _COM_Outptr_  void **ppvObject);
        
        ULONG ( STDMETHODCALLTYPE *AddRef )( 
            IEngine * This);
        
        ULONG ( STDMETHODCALLTYPE *Release )( 
            IEngine * This);
        
        HRESULT ( STDMETHODCALLTYPE *SpeedUp )( 
            IEngine * This);
        
        HRESULT ( STDMETHODCALLTYPE *GetMaxSpeed )( 
            IEngine * This,
            /* [out] */ int *maxSpeed);
        
        HRESULT ( STDMETHODCALLTYPE *GetCurSpeed )( 
            IEngine * This,
            /* [out] */ int *curSpeed);
        
        END_INTERFACE
    } IEngineVtbl;

    interface IEngine
    {
        CONST_VTBL struct IEngineVtbl *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IEngine_QueryInterface(This,riid,ppvObject)	\
    ( (This)->lpVtbl -> QueryInterface(This,riid,ppvObject) ) 

#define IEngine_AddRef(This)	\
    ( (This)->lpVtbl -> AddRef(This) ) 

#define IEngine_Release(This)	\
    ( (This)->lpVtbl -> Release(This) ) 


#define IEngine_SpeedUp(This)	\
    ( (This)->lpVtbl -> SpeedUp(This) ) 

#define IEngine_GetMaxSpeed(This,maxSpeed)	\
    ( (This)->lpVtbl -> GetMaxSpeed(This,maxSpeed) ) 

#define IEngine_GetCurSpeed(This,curSpeed)	\
    ( (This)->lpVtbl -> GetCurSpeed(This,curSpeed) ) 

#endif /* COBJMACROS */


#endif 	/* C style interface */




#endif 	/* __IEngine_INTERFACE_DEFINED__ */


#ifndef __IStats_INTERFACE_DEFINED__
#define __IStats_INTERFACE_DEFINED__

/* interface IStats */
/* [unique][helpstring][dual][uuid][object] */ 


EXTERN_C const IID IID_IStats;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("46073856-735C-4D0A-8C30-80BF4ECFC2DC")
    IStats : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE DisplayStats( void) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE GetPetName( 
            /* [out] */ BSTR *petName) = 0;
        
    };
    
    
#else 	/* C style interface */

    typedef struct IStatsVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE *QueryInterface )( 
            IStats * This,
            /* [in] */ REFIID riid,
            /* [annotation][iid_is][out] */ 
            _COM_Outptr_  void **ppvObject);
        
        ULONG ( STDMETHODCALLTYPE *AddRef )( 
            IStats * This);
        
        ULONG ( STDMETHODCALLTYPE *Release )( 
            IStats * This);
        
        HRESULT ( STDMETHODCALLTYPE *DisplayStats )( 
            IStats * This);
        
        HRESULT ( STDMETHODCALLTYPE *GetPetName )( 
            IStats * This,
            /* [out] */ BSTR *petName);
        
        END_INTERFACE
    } IStatsVtbl;

    interface IStats
    {
        CONST_VTBL struct IStatsVtbl *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IStats_QueryInterface(This,riid,ppvObject)	\
    ( (This)->lpVtbl -> QueryInterface(This,riid,ppvObject) ) 

#define IStats_AddRef(This)	\
    ( (This)->lpVtbl -> AddRef(This) ) 

#define IStats_Release(This)	\
    ( (This)->lpVtbl -> Release(This) ) 


#define IStats_DisplayStats(This)	\
    ( (This)->lpVtbl -> DisplayStats(This) ) 

#define IStats_GetPetName(This,petName)	\
    ( (This)->lpVtbl -> GetPetName(This,petName) ) 

#endif /* COBJMACROS */


#endif 	/* C style interface */




#endif 	/* __IStats_INTERFACE_DEFINED__ */


#ifndef __ICreateCar_INTERFACE_DEFINED__
#define __ICreateCar_INTERFACE_DEFINED__

/* interface ICreateCar */
/* [unique][helpstring][dual][uuid][object] */ 


EXTERN_C const IID IID_ICreateCar;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("F426AD9B-76F0-46BE-95DE-0B22C84EB89D")
    ICreateCar : public IUnknown
    {
    public:
        virtual HRESULT STDMETHODCALLTYPE SetPetName( 
            /* [in] */ BSTR maxSpeed) = 0;
        
        virtual HRESULT STDMETHODCALLTYPE SetMaxSpeed( 
            /* [in] */ int curSpeed) = 0;
        
    };
    
    
#else 	/* C style interface */

    typedef struct ICreateCarVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE *QueryInterface )( 
            ICreateCar * This,
            /* [in] */ REFIID riid,
            /* [annotation][iid_is][out] */ 
            _COM_Outptr_  void **ppvObject);
        
        ULONG ( STDMETHODCALLTYPE *AddRef )( 
            ICreateCar * This);
        
        ULONG ( STDMETHODCALLTYPE *Release )( 
            ICreateCar * This);
        
        HRESULT ( STDMETHODCALLTYPE *SetPetName )( 
            ICreateCar * This,
            /* [in] */ BSTR maxSpeed);
        
        HRESULT ( STDMETHODCALLTYPE *SetMaxSpeed )( 
            ICreateCar * This,
            /* [in] */ int curSpeed);
        
        END_INTERFACE
    } ICreateCarVtbl;

    interface ICreateCar
    {
        CONST_VTBL struct ICreateCarVtbl *lpVtbl;
    };

    

#ifdef COBJMACROS


#define ICreateCar_QueryInterface(This,riid,ppvObject)	\
    ( (This)->lpVtbl -> QueryInterface(This,riid,ppvObject) ) 

#define ICreateCar_AddRef(This)	\
    ( (This)->lpVtbl -> AddRef(This) ) 

#define ICreateCar_Release(This)	\
    ( (This)->lpVtbl -> Release(This) ) 


#define ICreateCar_SetPetName(This,maxSpeed)	\
    ( (This)->lpVtbl -> SetPetName(This,maxSpeed) ) 

#define ICreateCar_SetMaxSpeed(This,curSpeed)	\
    ( (This)->lpVtbl -> SetMaxSpeed(This,curSpeed) ) 

#endif /* COBJMACROS */


#endif 	/* C style interface */




#endif 	/* __ICreateCar_INTERFACE_DEFINED__ */


#ifndef __ICar_INTERFACE_DEFINED__
#define __ICar_INTERFACE_DEFINED__

/* interface ICar */
/* [unique][nonextensible][dual][uuid][object] */ 


EXTERN_C const IID IID_ICar;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("a3194885-8306-4df2-ab44-abcecec6c700")
    ICar : public IDispatch
    {
    public:
    };
    
    
#else 	/* C style interface */

    typedef struct ICarVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE *QueryInterface )( 
            ICar * This,
            /* [in] */ REFIID riid,
            /* [annotation][iid_is][out] */ 
            _COM_Outptr_  void **ppvObject);
        
        ULONG ( STDMETHODCALLTYPE *AddRef )( 
            ICar * This);
        
        ULONG ( STDMETHODCALLTYPE *Release )( 
            ICar * This);
        
        HRESULT ( STDMETHODCALLTYPE *GetTypeInfoCount )( 
            ICar * This,
            /* [out] */ UINT *pctinfo);
        
        HRESULT ( STDMETHODCALLTYPE *GetTypeInfo )( 
            ICar * This,
            /* [in] */ UINT iTInfo,
            /* [in] */ LCID lcid,
            /* [out] */ ITypeInfo **ppTInfo);
        
        HRESULT ( STDMETHODCALLTYPE *GetIDsOfNames )( 
            ICar * This,
            /* [in] */ REFIID riid,
            /* [size_is][in] */ LPOLESTR *rgszNames,
            /* [range][in] */ UINT cNames,
            /* [in] */ LCID lcid,
            /* [size_is][out] */ DISPID *rgDispId);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE *Invoke )( 
            ICar * This,
            /* [annotation][in] */ 
            _In_  DISPID dispIdMember,
            /* [annotation][in] */ 
            _In_  REFIID riid,
            /* [annotation][in] */ 
            _In_  LCID lcid,
            /* [annotation][in] */ 
            _In_  WORD wFlags,
            /* [annotation][out][in] */ 
            _In_  DISPPARAMS *pDispParams,
            /* [annotation][out] */ 
            _Out_opt_  VARIANT *pVarResult,
            /* [annotation][out] */ 
            _Out_opt_  EXCEPINFO *pExcepInfo,
            /* [annotation][out] */ 
            _Out_opt_  UINT *puArgErr);
        
        END_INTERFACE
    } ICarVtbl;

    interface ICar
    {
        CONST_VTBL struct ICarVtbl *lpVtbl;
    };

    

#ifdef COBJMACROS


#define ICar_QueryInterface(This,riid,ppvObject)	\
    ( (This)->lpVtbl -> QueryInterface(This,riid,ppvObject) ) 

#define ICar_AddRef(This)	\
    ( (This)->lpVtbl -> AddRef(This) ) 

#define ICar_Release(This)	\
    ( (This)->lpVtbl -> Release(This) ) 


#define ICar_GetTypeInfoCount(This,pctinfo)	\
    ( (This)->lpVtbl -> GetTypeInfoCount(This,pctinfo) ) 

#define ICar_GetTypeInfo(This,iTInfo,lcid,ppTInfo)	\
    ( (This)->lpVtbl -> GetTypeInfo(This,iTInfo,lcid,ppTInfo) ) 

#define ICar_GetIDsOfNames(This,riid,rgszNames,cNames,lcid,rgDispId)	\
    ( (This)->lpVtbl -> GetIDsOfNames(This,riid,rgszNames,cNames,lcid,rgDispId) ) 

#define ICar_Invoke(This,dispIdMember,riid,lcid,wFlags,pDispParams,pVarResult,pExcepInfo,puArgErr)	\
    ( (This)->lpVtbl -> Invoke(This,dispIdMember,riid,lcid,wFlags,pDispParams,pVarResult,pExcepInfo,puArgErr) ) 


#endif /* COBJMACROS */


#endif 	/* C style interface */




#endif 	/* __ICar_INTERFACE_DEFINED__ */



#ifndef __CarComponentLib_LIBRARY_DEFINED__
#define __CarComponentLib_LIBRARY_DEFINED__

/* library CarComponentLib */
/* [version][uuid] */ 


EXTERN_C const IID LIBID_CarComponentLib;

EXTERN_C const CLSID CLSID_Car;

#ifdef __cplusplus

class DECLSPEC_UUID("c779aa99-9945-40b6-b144-26b1cfe0a12c")
Car;
#endif
#endif /* __CarComponentLib_LIBRARY_DEFINED__ */

/* Additional Prototypes for ALL interfaces */

unsigned long             __RPC_USER  BSTR_UserSize(     unsigned long *, unsigned long            , BSTR * ); 
unsigned char * __RPC_USER  BSTR_UserMarshal(  unsigned long *, unsigned char *, BSTR * ); 
unsigned char * __RPC_USER  BSTR_UserUnmarshal(unsigned long *, unsigned char *, BSTR * ); 
void                      __RPC_USER  BSTR_UserFree(     unsigned long *, BSTR * ); 

unsigned long             __RPC_USER  BSTR_UserSize64(     unsigned long *, unsigned long            , BSTR * ); 
unsigned char * __RPC_USER  BSTR_UserMarshal64(  unsigned long *, unsigned char *, BSTR * ); 
unsigned char * __RPC_USER  BSTR_UserUnmarshal64(unsigned long *, unsigned char *, BSTR * ); 
void                      __RPC_USER  BSTR_UserFree64(     unsigned long *, BSTR * ); 

/* end of Additional Prototypes */

#ifdef __cplusplus
}
#endif

#endif


