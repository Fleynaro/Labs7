using System;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Runtime.InteropServices.ComTypes;

namespace CarCsClient
{
    internal static class Program
    {
        [STAThread]
        private static void Main()
        {
            var carType = Type.GetTypeFromProgID("CarComponent.Car");
            dynamic car = Activator.CreateInstance(carType!);

            //var interfaces = carType.GetInterfaces();
            //car.DisplayStats();
            //carType.InvokeMember("SpeedUp", BindingFlags.InvokeMethod, null, car, new object[] { });

            var iDisp = (IDispatch) car;
            var count = iDisp!.GetTypeInfoCount();
            Console.WriteLine($"Type info count is {count}.");

            //var dispId = new int []{};
            //iDisp.GetIDsOfNames(Guid.Empty, new [] {"SetMaxSpeed"}, 1, 0x800, dispId);

            var iEngine = (IEngine) car;
            iEngine!.SpeedUp();
            var speed = 0;
            iEngine.GetCurSpeed(ref speed);
            Console.WriteLine($"Car speed is {speed}.");

            var iCrCar = (ICreateCar) car;
            Console.Write("Enter pet name: ");
            iCrCar!.SetPetName(Console.ReadLine());

            var iStats = (IStats) car;
            string petName = null;
            iStats!.GetPetName(ref petName);
            Console.WriteLine($"Pet name is '{petName}'.");
        }
    }

    [ComImport, InterfaceType(ComInterfaceType.InterfaceIsIUnknown), Guid("F426AD9B-76F0-46BE-95DE-0B22C84EB89D")]
    public interface ICreateCar
    {
        void SetPetName(string petName);
        void SetMaxSpeed(int maxSp);
    }

    [ComImport, InterfaceType(ComInterfaceType.InterfaceIsIUnknown), Guid("46073856-735C-4D0A-8C30-80BF4ECFC2DC")]
    public interface IStats
    {
        void DisplayStats();
        void GetPetName(ref string petName);
    }

    [ComImport, InterfaceType(ComInterfaceType.InterfaceIsIUnknown), Guid("9E37E773-86BD-4DD5-8F74-30536AAFC3B3")]
    public interface IEngine
    {
        void SpeedUp();
        void GetMaxSpeed(ref int curSpeed);
        void GetCurSpeed(ref int maxSpeed);
    }

    [ComImport, InterfaceType(ComInterfaceType.InterfaceIsIUnknown), Guid("00020400-0000-0000-C000-000000000046")]
    public interface IDispatch
    {
        int GetTypeInfoCount();

        [return: MarshalAs(UnmanagedType.Interface)]
        ITypeInfo GetTypeInfo([In, MarshalAs(UnmanagedType.U4)] int iTInfo, [In, MarshalAs(UnmanagedType.U4)] int lcid);

        void GetIDsOfNames([In] ref Guid riid, [In, MarshalAs(UnmanagedType.LPArray)] string[] rgszNames,
            [In, MarshalAs(UnmanagedType.U4)] int cNames, [In, MarshalAs(UnmanagedType.U4)] int lcid,
            [Out, MarshalAs(UnmanagedType.LPArray)]
            int[] rgDispId);
    }
}
