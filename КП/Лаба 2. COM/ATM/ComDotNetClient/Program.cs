using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace ComDotNetClient
{
    class Program
    {
        [ComVisible(true)]
        [ComImport, Guid("B18544E2-DFCE-4A29-8EAE-E13F4A5FB575")]
        public class ATM
        {
        }

        [ComImport, InterfaceType(ComInterfaceType.InterfaceIsIUnknown), Guid("D69D0D10-76CA-40D0-A6AA-87C04B0C0662")]
        public interface ICreateAtm
        {
            void SetLocationAddr(string addr);
            void SetMoney(int money);
        }

        [ComImport, InterfaceType(ComInterfaceType.InterfaceIsIUnknown), Guid("56789D0E-242E-459F-9334-EA8BCBE14DE8")]
        public interface IStats
        {
            void DisplayStats();
            void GetLocationAddr(ref string addr);
        }

        [ComImport, InterfaceType(ComInterfaceType.InterfaceIsIUnknown), Guid("A2D3FDF3-7083-4A2F-86B5-BC238CEF49EF")]
        public interface IBalance
        {
            void PayMoney(int accountId, int money);
            void WithdrawMoney(int accountId, int money);
            void GetBalance(int accountId, ref int balance);
        }

        [STAThreadAttribute]
        static void Main(string[] args)
        {
            ATM atm = new ATM();

            ICreateAtm iCrAtn = (ICreateAtm)atm;
            Console.WriteLine("Напишите адрес банкомата: ");
            iCrAtn.SetLocationAddr(Console.ReadLine());

            Console.ReadKey();
        }
    }
}
