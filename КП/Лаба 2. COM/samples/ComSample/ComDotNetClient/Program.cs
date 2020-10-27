using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace ComDotNetClient
{
    [ComVisible(true)]
    [ComImport, Guid("7AD2D539-EE35-11d2-B8DE-0020781238D4")]
    public class Car
    {
    }

    //{960EFB55-5A90-45F8-8C6C-12849DCF2AC1}
    //DEFINE_GUID(IID_ICreateCar,
    //0x960efb55, 0x5a90, 0x45f8, 0x8c, 0x6c, 0x12, 0x84, 0x9d, 0xcf, 0x2a, 0xc1);
    [ComImport, InterfaceType(ComInterfaceType.InterfaceIsIUnknown), Guid("960EFB55-5A90-45F8-8C6C-12849DCF2AC1")]
    public interface ICreateCar
    {
        void SetPetName(string petName);
        void SetMaxSpeed(int maxSp);
    }

    // {4601678B-6D77-446F-8782-CC054440F81F}
    //DEFINE_GUID(IID_IStats,
    //0x4601678b, 0x6d77, 0x446f, 0x87, 0x82, 0xcc, 0x5, 0x44, 0x40, 0xf8, 0x1f);
    [ComImport, InterfaceType(ComInterfaceType.InterfaceIsIUnknown), Guid("4601678B-6D77-446F-8782-CC054440F81F")]
    public interface IStats
    {
        void DisplayStats();
        void GetPetName(ref string petName);
    }
    [ComImport, InterfaceType(ComInterfaceType.InterfaceIsIUnknown), Guid("0010D07E-8AF4-4436-A4D1-A2F57C6A4DA4")]
    public interface IEngine
    {
        void SpeedUp();
        void GetMaxSpeed(ref int curSpeed);
        void GetCurSpeed(ref int maxSpeed);
    }


    class Program
    {
        static void Main(string[] args)
        {
            Car myCar = new Car();
            
            ICreateCar iCrCar = (ICreateCar)myCar;
            Console.WriteLine("Напишите имя: ");
            iCrCar.SetPetName(Console.ReadLine());
        }
    }
}
