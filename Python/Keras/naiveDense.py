import tensorflow as tf

class NaiveDense:
    def __init__(self, input_size, output_size, activation):
        self.activation = activation

        w_initial_value = tf.random.uniform((input_size, output_size), minval=0, maxval=0.1)
        self.W = tf.Variable(w_initial_value)

        b_initial_value = tf.zeros((output_size, ))
        self.b = tf.Variable(b_initial_value)

    def __call__(self, inputs):
        # print("inputs.shape",inputs.shape)
        # print("self.W.shape", self.W.shape)
        # exit()
        
        return self.activation(tf.matmul(inputs, self.W)+self.b)

    @property
    def weights(self):
        return [self.W, self.b]
