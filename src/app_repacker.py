import os

from utils.filesystem_utils import get_bucket_filesystem


def main() -> None:
    bucket_fs = get_bucket_filesystem()
    gnome_bucket = os.environ.get("GNOME_BUCKET", "3dgnome-landing-zone")
    model_repository_bucket = os.environ.get("MODEL_REPOSITORY_BUCKET", "model-repository")




if __name__ == '__main__':
    main()
