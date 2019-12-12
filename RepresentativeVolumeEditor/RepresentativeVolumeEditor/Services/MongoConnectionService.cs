using MongoDB.Driver;
using System;
using System.Collections.Generic;
using System.Text;

namespace RepresentativeVolumeEditor.Services
{
    static class MongoConnectionService
    {
        public static MongoClient Connect(string connectionString)
        {
            if (string.IsNullOrEmpty(connectionString))
            {
                throw new ArgumentNullException(nameof(connectionString));
            }

            var dbClient = new MongoClient(connectionString);

            return dbClient;
        }
    }
}
