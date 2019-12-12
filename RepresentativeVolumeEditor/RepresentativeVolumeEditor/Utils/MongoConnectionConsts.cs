using System;
using System.Collections.Generic;
using System.Text;

namespace RepresentativeVolumeEditor.Utils
{
    static class MongoConnectionConsts
    {
        public static readonly string ConnectionString = "mongodb+srv://m001-student:m001-mongodb-basics@sandbox-0tzkj.mongodb.net/test?retryWrites=true&w=majority";
        public static readonly string DatabaseName = "thermalConductivityIdentification";
        public static readonly string CollectionName = "representativeVolumeMaterials";
    }
}
