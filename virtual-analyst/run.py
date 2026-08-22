"""Starts the development server.

    py run.py

It is not a production server - and does not need to be: the app runs on the
analyst's machine, with one user. `threaded=True` matters because the screen
polls while the background worker is busy.
"""

from app import config
from app import create_app

app = create_app()

if __name__ == "__main__":
    app.run(
        host="127.0.0.1",
        port=5000,
        debug=config.DEBUG,
        threaded=True,
        # Flask's reloader starts TWO processes - and the background worker in
        # the second one would be a ghost worker draining the same queue.
        use_reloader=False,
    )
