model_cov = keras.models.Sequential([
      keras.layers.Conv2D(64,7,activation="relu",padding="same",input_shape=[64,64,3]),
      keras.layers.MaxPooling2D(2),
      keras.layers.Conv2D(128,3,activation="relu",padding="same"),
      keras.layers.Conv2D(128,3,activation="relu",padding="same"),
      keras.layers.MaxPooling2D(2),
      keras.layers.Conv2D(256,3,activation="relu",padding="same"),
      keras.layers.Conv2D(256,3,activation="relu",padding="same"),
      keras.layers.MaxPooling2D(2),
      keras.layers.Flatten(),
      layers.Dense(128, kernel_initializer="he_normal", name="layer1",kernel_regularizer=tf.keras.regularizers.L2(0.01)),
      keras.layers.LeakyReLU(alpha=0.2),
      layers.Dense(64, kernel_initializer="he_normal", name="layer2",kernel_regularizer=tf.keras.regularizers.L2(0.01)),
      keras.layers.LeakyReLU(alpha=0.2),
      keras.layers.Dense(10,activation="softmax")
])

