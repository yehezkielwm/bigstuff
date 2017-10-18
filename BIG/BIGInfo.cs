using System;
using System.Drawing;
using Grasshopper.Kernel;

namespace BIG
{
    public class BIGInfo : GH_AssemblyInfo
    {
        public override string Name
        {
            get
            {
                return "BIG";
            }
        }
        public override Bitmap Icon
        {
            get
            {
                //Return a 24x24 pixel bitmap to represent this GHA library.
                return null;
            }
        }
        public override string Description
        {
            get
            {
                //Return a short string describing the purpose of this GHA library.
                return "Useful tools for BIG Ideas projects and analysis";
            }
        }
        public override Guid Id
        {
            get
            {
                return new Guid("24025dd1-7598-4e08-9697-bea6d4acbb63");
            }
        }

        public override string AuthorName
        {
            get
            {
                //Return a string identifying you or your company.
                return "BIG";
            }
        }
        public override string AuthorContact
        {
            get
            {
                //Return a string representing your preferred contact details.
                return "yewi@big.dk, krne@big.dk, tore@big.dk";
            }
        }
    }
}
