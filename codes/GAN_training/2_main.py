import os

from WGAN_GP import WGAN_GP

# from utils import show_all_variables
from utils import check_folder


# import tensorflow as tf
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()

import argparse

"""parsing and configuration"""
def parse_args():
    desc = "Tensorflow implementation of GAN collections"
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('--gan_type', type=str, default='WGAN_GP',
                        choices=['GAN', 'CGAN', 'infoGAN', 'ACGAN', 'EBGAN', 'BEGAN', 'WGAN', 'WGAN_GP', 'DRAGAN', 'LSGAN', 'VAE', 'CVAE'],
                        help='The type of GAN') # required=True)
    parser.add_argument('--dataset', type=str, default='BLA_10k', choices=['PFC', 'CPU', 'HIP', 'NAC', 'VTA'],
                        help='The name of dataset')
    parser.add_argument('--epoch', type=int, default=50001, help='The number of epochs to run')
    # Confusing, in the body of paper, author said batch_size was 32. but in supplementary, Its 30
    parser.add_argument('--batch_size', type=int, default=32, help='The size of batch')
    parser.add_argument('--z_dim', type=int, default=100, help='Dimension of noise vector')
    parser.add_argument('--checkpoint_dir', type=str, default='checkpoint',
                        help='Directory name to save the checkpoints')
    parser.add_argument('--result_dir', type=str, default='results',
                        help='Directory name to save the generated images')
    parser.add_argument('--log_dir', type=str, default='logs',
                        help='Directory name to save training logs')

    # some arguments added from original code

    parser.add_argument('--inflate_to_size', type=int, default=800, help='inflate_to_size in Generator')
    parser.add_argument('--gex_size', type=int, default=13557, help='gex_size; gene expression size of data')
    parser.add_argument('--disc_internal_size', type=int, default=400, help='disc_internal_size in Discriminator')
    parser.add_argument('--gpu_id', type=int, default=0, help='using gpu_id')


    return check_args(parser.parse_args())

"""checking arguments"""
def check_args(args):
    # --checkpoint_dir
    check_folder(args.checkpoint_dir)

    # --result_dir
    check_folder(args.result_dir)

    # --result_dir
    check_folder(args.log_dir)

    # --epoch
    assert args.epoch >= 1, 'number of epochs must be larger than or equal to one'

    # --batch_size
    assert args.batch_size >= 1, 'batch size must be larger than or equal to one'

    # --z_dim
    assert args.z_dim >= 1, 'dimension of noise vector must be larger than or equal to one'

    return args

"""main"""
def main():
    # parse arguments
    args = parse_args()
    if args is None:
      exit()

    # open session
    # in case of other Models
    # models = [GAN, CGAN, infoGAN, ACGAN, EBGAN, WGAN, WGAN_GP, DRAGAN,
    #           LSGAN, BEGAN, VAE, CVAE]
    models = [WGAN_GP]
    os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu_id)

    with tf.Session(config=tf.ConfigProto(allow_soft_placement=True)) as sess:

        # declare instance for GAN

        gan = None
        for model in models:
            if args.gan_type == model.model_name:
                gan = model(sess,
                            epoch=args.epoch,
                            batch_size=args.batch_size,
                            z_dim=args.z_dim,
                            dataset_name=args.dataset,
                            checkpoint_dir=args.checkpoint_dir,
                            result_dir=args.result_dir,
                            log_dir=args.log_dir,
                            # added arguments
                            inflate_to_size=args.inflate_to_size,
                            gex_size=args.gex_size,
                            disc_internal_size=args.disc_internal_size
			    )
        if gan is None:
            raise Exception("[!] There is no option for " + args.gan_type)

        # build graph
        gan.build_model()

        # show network architecture
        # show_all_variables()

        # launch the graph in a session
        gan.train()
        print(" [*] Training finished!")

        # visualize learned generator
        gan.visualize_results(args.epoch-1)
        print(" [*] Testing finished!")

if __name__ == '__main__':
    main()
