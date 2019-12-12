using MongoDB.Bson;
using MongoDB.Driver;
using RepresentativeVolumeEditor.Utils;
using System;
using System.Collections.Generic;
using System.Text;

namespace RepresentativeVolumeEditor.Services
{
    static class RepresentativeVolumeMaterialsService
    {
        private static IMongoCollection<BsonDocument> collection;
        
        public static List<Material> GetMaterials()
        {
            if (collection == null)
            {
                LoadCollection();
            }

            var result = new List<Material>();
            var rawMaterials = collection.Find<BsonDocument>(new BsonDocument()).ToList();
            foreach (var item in rawMaterials)
            {
                var material = new Material(item);
                result.Add(material);
            }

            return result;
        }

        public static void AddMaterial(Material material)
        {
            if (material == null)
            {
                throw new ArgumentNullException(nameof(material));
            }

            if (collection == null)
            {
                LoadCollection();
            }

            var bsonMaterial = material.ToBsonDocument();
            collection.InsertOne(bsonMaterial);
        }

        private static void LoadCollection()
        {
            var dbClient = MongoConnectionService.Connect(MongoConnectionConsts.ConnectionString);
            if (dbClient == null)
            {
                throw new NullReferenceException(nameof(dbClient));
            }

            var database = dbClient.GetDatabase(MongoConnectionConsts.DatabaseName);
            if (database == null)
            {
                throw new NullReferenceException(nameof(database));
            }

            collection = database.GetCollection<BsonDocument>(MongoConnectionConsts.CollectionName);
            if (collection == null)
            {
                throw new NullReferenceException(nameof(collection));
            }
        }
    }
}
