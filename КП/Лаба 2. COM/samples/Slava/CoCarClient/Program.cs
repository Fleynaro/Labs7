using System;
using System.Runtime.InteropServices;

namespace CoCarClientCs {
    class Program {

        [ComVisible(true)]
        [ComImport, Guid("52060B95-5EDF-46B4-A487-C153FEA37D79")]
        public class Car {
        }

        [ComImport, InterfaceType(ComInterfaceType.InterfaceIsIUnknown), Guid("272A1A41-442E-4863-98F6-49815DF857ED")]
        public interface ICreateCar {
            void SetPetName(string petName);
            void SetMaxSpeed(int maxSp);
        }

        [ComImport, InterfaceType(ComInterfaceType.InterfaceIsIUnknown), Guid("F57C14BF-FC1B-44E5-98EF-2A3F55F63F13")]
        public interface IStats {
            void DisplayStats();
            void GetPetName(ref string petName);
        }

        [ComImport, InterfaceType(ComInterfaceType.InterfaceIsIUnknown), Guid("796E603B-1FB0-4A6E-A3C0-67FD20163E3B")]
        public interface IEngine {
            void SpeedUp();
            void GetMaxSpeed(ref int curSpeed);
            void GetCurSpeed(ref int maxSpeed);
        }


        [STAThreadAttribute]
        static void Main(string[] args) {
            Car myCar = new Car();
            ICreateCar iCrCar = (ICreateCar)myCar;
            Console.WriteLine("Напишите имя: ");
            iCrCar.SetPetName(Console.ReadLine());
            iCrCar.SetPetName(Console.ReadLine());
        }
    }
}

