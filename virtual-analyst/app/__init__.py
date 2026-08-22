"""Virtual Analyst - app factory.

No access control, by design: this is a bench tool, running on the machine of
whoever does the analysis. Do not expose it on a network without putting
authentication in front of it.
"""

from __future__ import annotations

from flask import Flask, render_template

from app import config


def create_app() -> Flask:
    app = Flask(__name__)
    app.config.update(
        SECRET_KEY=config.SECRET_KEY,
        MAX_CONTENT_LENGTH=config.MAX_UPLOAD_MB * 1024 * 1024,
        # Without this, editing a template means restarting the server.
        TEMPLATES_AUTO_RELOAD=config.DEBUG,
        JSON_SORT_KEYS=False,
    )

    from app.blueprints.method import bp as method_bp
    from app.blueprints.results import bp as results_bp
    from app.blueprints.samples import bp as samples_bp
    from app.blueprints.tasks import bp as tasks_bp

    app.register_blueprint(samples_bp)
    app.register_blueprint(method_bp)
    app.register_blueprint(tasks_bp)
    app.register_blueprint(results_bp)

    @app.route("/")
    def home():
        return render_template("home.html")

    @app.errorhandler(413)
    def upload_too_large(_error):
        return (
            render_template(
                "error.html",
                title="Upload too large",
                message=(
                    f"The current ceiling is {config.MAX_UPLOAD_MB} MB for the whole request "
                    f"(the sum of every file sent at once). Upload fewer files per batch or "
                    f"raise MAX_UPLOAD_MB in the .env."
                ),
            ),
            413,
        )

    return app
